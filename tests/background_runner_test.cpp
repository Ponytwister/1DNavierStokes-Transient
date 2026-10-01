#include <background_runner.h>
#include <gtest/gtest.h>
#include <future>

namespace {
using namespace tsensor_workflow;

// Synchronize at a real progress boundary, without timing-dependent sleeps.
struct gate {
    std::promise<void> entered;
    std::promise<void> released;
    std::shared_future<void> release = released.get_future().share();
    bool opened = false;
    void open() { if (!opened) { opened = true; released.set_value(); } }
    ~gate() { open(); }
    progress_callback callback() {
        return [this](const progress_event& event) {
            if (event.kind == event_kind::started && event.action == operation::open_database) {
                entered.set_value();
                release.wait();
            }
        };
    }
};

void expect_cancelled(const std::exception_ptr& failure)
{
    ASSERT_TRUE(failure);
    try { std::rethrow_exception(failure); }
    catch (const workflow_error& error) { EXPECT_EQ(error.code, error_code::cancelled); }
}

TEST(BackgroundRunner, SingleActiveRunCancellationAndReuse)
{
    background_runner runner;
    gate blocked; // Releases before runner destruction, including assertion failures.
    auto entered = blocked.entered.get_future();
    runner.start("unused-cancelled-test.db", blocked.callback());
    ASSERT_EQ(entered.wait_for(std::chrono::seconds(5)), std::future_status::ready);
    EXPECT_EQ(runner.status(), background_state::running);
    EXPECT_THROW(runner.start("other.db"), workflow_error);
    EXPECT_THROW(runner.take_result(), workflow_error);
    auto events = runner.drain_events();
    ASSERT_EQ(events.size(), 1u);
    EXPECT_EQ(events.front().action, operation::open_database);
    EXPECT_TRUE(runner.drain_events().empty());
    EXPECT_TRUE(runner.request_cancel());
    EXPECT_FALSE(runner.request_cancel());
    blocked.open();
    runner.wait();
    EXPECT_EQ(runner.status(), background_state::cancelled);
    EXPECT_THROW(runner.start("other.db"), workflow_error);
    auto outcome = runner.take_result();
    EXPECT_FALSE(outcome.session);
    EXPECT_FALSE(outcome.result);
    expect_cancelled(outcome.failure);
    EXPECT_EQ(runner.status(), background_state::idle);
    EXPECT_FALSE(runner.request_cancel());
    EXPECT_THROW(runner.take_result(), workflow_error);

    gate again;
    auto second = again.entered.get_future();
    runner.start("unused-cancelled-test.db", again.callback());
    ASSERT_EQ(second.wait_for(std::chrono::seconds(5)), std::future_status::ready);
    EXPECT_TRUE(runner.request_cancel());
    again.open();
    runner.wait();
    expect_cancelled(runner.take_result().failure);
}

TEST(BackgroundRunner, MissingDatabaseFailureIsRetained)
{
    background_runner runner;
    runner.start(std::filesystem::path(TSENSOR_FIXTURE_DIR) / "missing-directory" / "absent.db");
    runner.wait();
    EXPECT_EQ(runner.status(), background_state::failed);
    EXPECT_FALSE(runner.request_cancel());
    const auto events = runner.drain_events();
    ASSERT_GE(events.size(), 2u);
    EXPECT_EQ(events.back().error, error_code::database);
    auto outcome = runner.take_result();
    ASSERT_TRUE(outcome.failure);
    EXPECT_FALSE(outcome.session);
    EXPECT_FALSE(outcome.result);
    try { std::rethrow_exception(outcome.failure); }
    catch (const workflow_error& error) {
        EXPECT_EQ(error.code, error_code::database);
        EXPECT_EQ(error.action, operation::open_database);
    }
}

TEST(BackgroundRunner, DestructionJoinsWorker)
{
    auto runner = std::make_unique<background_runner>();
    gate blocked;
    auto entered = blocked.entered.get_future();
    runner->start("unused-cancelled-test.db", blocked.callback());
    ASSERT_EQ(entered.wait_for(std::chrono::seconds(5)), std::future_status::ready);
    // Release the callback concurrently with destruction; the destructor must
    // not return with a worker still accessing runner members.
    std::jthread release([&] { blocked.open(); });
    runner.reset();
    release.join();
}

TEST(Cancellation, CheckpointsPrecedeDatabaseAndModelAccess)
{
    std::stop_source stop;
    parameters_t p;
    p.cancellation = stop.get_token();
    p.active_operation = operation::solve;
    stop.request_stop();
    alglib::real_1d_array controls, residuals;
    try { alglib_solver(controls, residuals, &p); FAIL(); }
    catch (const workflow_error& error) { EXPECT_EQ(error.code, error_code::cancelled); }
    try { read_model_parameters_from_db(p, nullptr); FAIL(); }
    catch (const workflow_error& error) { EXPECT_EQ(error.code, error_code::cancelled); }
    EXPECT_EQ(p.iterations, 0);
}

TEST(Cancellation, SessionLoadRejectsStoppedTokenWithoutDatabaseWrites)
{
    std::stop_source stop;
    run_session session(":memory:", {}, stop.get_token());
    stop.request_stop();
    try { session.load_inputs(); FAIL(); }
    catch (const workflow_error& error) {
        EXPECT_EQ(error.code, error_code::cancelled);
        EXPECT_EQ(error.action, operation::load_inputs);
    }
    EXPECT_EQ(sqlite3_total_changes(session.database()), 0);
    EXPECT_THROW(session.export_results("unused-output"), workflow_error);
}

TEST(Cancellation, AlglibUnwindsCancelledCallbackAndCanRunAgain)
{
    // Independent scalar least-squares case: residual x-2 has minimum x=2.
    for (bool cancel : {true, false}) {
        parameters_t p;
        std::stop_source stop;
        p.cancellation = stop.get_token();
        p.active_operation = operation::solve;
        if (cancel) stop.request_stop();
        alglib::real_1d_array x = "[0]";
        alglib::minlmstate state;
        alglib::minlmreport report;
        alglib::minlmcreatev(1, 1, x, 0.0001, state);
        alglib::minlmsetcond(state, 1e-12, 100);
        const auto residual = [](const alglib::real_1d_array& x, alglib::real_1d_array& fi, void* ptr) {
            check_cancellation(*static_cast<parameters_t*>(ptr));
            fi[0] = x[0] - 2.0;
        };
        if (cancel) {
            try { alglib::minlmoptimize(state, +residual, nullptr, &p); FAIL(); }
            catch (const workflow_error& error) { EXPECT_EQ(error.code, error_code::cancelled); }
        } else {
            alglib::minlmoptimize(state, +residual, nullptr, &p);
            alglib::minlmresults(state, x, report);
            EXPECT_TRUE(std::isfinite(x[0]));
            EXPECT_NEAR(x[0], 2.0, 1e-12 + 1e-10 * 2.0);
        }
    }
}
} // namespace
