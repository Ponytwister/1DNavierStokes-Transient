#include <workflow.h>
#include <gtest/gtest.h>

#include <atomic>

namespace {
using namespace tsensor_workflow;

void execute(sqlite3* db, const char* sql)
{
    if (sqlite3_exec(db, sql, nullptr, nullptr, nullptr) != SQLITE_OK) {
        throw std::runtime_error(sqlite3_errmsg(db));
    }
}

double scalar(sqlite3* db, const char* sql)
{
    sqlite3_stmt* raw = nullptr;
    const int rc = sqlite3_prepare_v2(db, sql, -1, &raw, nullptr);
    std::unique_ptr<sqlite3_stmt, decltype(&sqlite3_finalize)> statement(raw, sqlite3_finalize);
    if (rc != SQLITE_OK || sqlite3_step(raw) != SQLITE_ROW) {
        throw std::runtime_error(sqlite3_errmsg(db));
    }
    return sqlite3_column_double(raw, 0);
}

class console_capture {
    bool active = true;
public:
    console_capture() {
        testing::internal::CaptureStdout();
        testing::internal::CaptureStderr();
    }
    std::pair<std::string, std::string> finish() {
        active = false;
        return {testing::internal::GetCapturedStdout(), testing::internal::GetCapturedStderr()};
    }
    ~console_capture() { if (active) { finish(); } }
};

TEST(Feedback, SessionReportsLifecycleAndStructuredDatabaseFailure)
{
    std::vector<progress_event> events;
    console_capture capture;
    run_session session(":memory:", [&](const progress_event& event) { events.push_back(event); });
    try {
        session.load_inputs(); // Deliberately absent schema.
        FAIL() << "Expected a database failure";
    } catch (const workflow_error& error) {
        EXPECT_EQ(error.action, operation::load_inputs);
        EXPECT_EQ(error.code, error_code::database);
        EXPECT_EQ(error.sqlite_code, SQLITE_ERROR);
    }
    const auto output = capture.finish();
    EXPECT_TRUE(output.first.empty());
    EXPECT_TRUE(output.second.empty());
    ASSERT_GE(events.size(), 4u);
    EXPECT_EQ(events[0].kind, event_kind::started);
    EXPECT_EQ(events[0].action, operation::open_database);
    EXPECT_EQ(events[1].kind, event_kind::completed);
    EXPECT_EQ(events[1].action, operation::open_database);
    EXPECT_EQ(events[2].kind, event_kind::started);
    EXPECT_EQ(events[2].action, operation::load_inputs);
    EXPECT_EQ(events.back().kind, event_kind::failed);
    EXPECT_EQ(events.back().error, error_code::database);
    EXPECT_EQ(events.back().sqlite_code, SQLITE_ERROR);
    EXPECT_EQ(session.state(), session_state::failed);
}

TEST(Feedback, NoObserverMeansSilentCoreAndSavingRemainsExplicit)
{
    run_result result;
    console_capture capture;
    {
        run_session session(":memory:");
        execute(session.database(), R"(
            CREATE TABLE model_controls (criterion TEXT, value TEXT);
            INSERT INTO model_controls VALUES ('universal_solve_for', 'width');
            CREATE TABLE alglib_input (
                VARIABLE TEXT, "INITIAL VALUE" REAL, "LOWER BOUND" REAL,
                "UPPER BOUND" REAL, SCALE REAL);
            INSERT INTO alglib_input VALUES ('width', 2.0, 0.1, 10.0, 1.0);
        )");
        auto& p = session.parameters();
        read_model_parameters_from_db(p, session.database());
        read_alglib_values_from_db(p, session.database());
        // Exercise the established optimizer-disabled branch, not a full solve.
        p.run_solver = false;
        p.total_window_size = 1;
        p.initial_values_alglib[0] = 3.0;
        result = run(p);
        EXPECT_DOUBLE_EQ(scalar(session.database(), "SELECT \"INITIAL VALUE\" FROM alglib_input;"), 2.0);
        save_fitted_parameters(p, session.database());
        EXPECT_DOUBLE_EQ(scalar(session.database(), "SELECT \"INITIAL VALUE\" FROM alglib_input;"), 3.0);
    }
    const auto output = capture.finish();
    EXPECT_TRUE(output.first.empty());
    EXPECT_TRUE(output.second.empty());
    EXPECT_FALSE(result.optimizer_ran);
    EXPECT_FALSE(result.termination_type.has_value());
    EXPECT_FALSE(result.optimizer_iterations.has_value());
    EXPECT_EQ(result.residual_evaluations, 0);
    ASSERT_EQ(result.parameters.size(), 1u);
    EXPECT_EQ(result.parameters[0].source, "global");
    EXPECT_EQ(result.parameters[0].name, "width");
    EXPECT_TRUE(std::isfinite(result.parameters[0].value));
    EXPECT_DOUBLE_EQ(result.parameters[0].value, 3.0); // Snapshot outlives its session.
}

TEST(Feedback, EvaluationCountsBelongToTheirOwnState)
{
    parameters_t first{}, second{};
    std::vector<progress_event> events;
    first.debug_level = second.debug_level = 3;
    first.active_operation = operation::solve;
    add_report(first, 3, "Solving");
    first.progress = [&](const progress_event& event) { events.push_back(event); };
    alglib::real_1d_array controls, residuals;
    // Empty experiment lists isolate callback bookkeeping from numerical solving.
    alglib_solver(controls, residuals, &first);
    alglib_solver(controls, residuals, &second);
    alglib_solver(controls, residuals, &first);
    ASSERT_EQ(events.size(), 2u);
    EXPECT_EQ(events[0].kind, event_kind::evaluation);
    EXPECT_EQ(events[0].action, operation::solve);
    EXPECT_EQ(events[0].evaluations, 1);
    EXPECT_EQ(events[1].evaluations, 2);
    EXPECT_EQ(second.iterations, 1);
    ASSERT_EQ(first.state.size(), 1u);
    EXPECT_EQ(first.state.back().debug_level, 3);
}

TEST(Feedback, ThrowingObserverIsDisconnectedWithoutChangingModelErrors)
{
    int calls = 0;
    run_session session(":memory:", [&](const progress_event&) {
        ++calls;
        throw std::runtime_error("display unavailable");
    });
    EXPECT_EQ(calls, 1);
    EXPECT_EQ(session.state(), session_state::empty);
    ASSERT_TRUE(session.progress_failure());
    EXPECT_THROW(std::rethrow_exception(session.progress_failure()), std::runtime_error);
    try {
        session.load_inputs();
        FAIL() << "Expected missing-schema error";
    } catch (const workflow_error& error) {
        EXPECT_EQ(error.code, error_code::database);
        EXPECT_EQ(error.action, operation::load_inputs);
    }
    EXPECT_EQ(calls, 1);
}

TEST(Feedback, WorkerDeliveriesAreSerialized)
{
    parameters_t p{};
    std::atomic<int> active{0}, overlaps{0};
    int delivered = 0;
    p.progress = [&](const progress_event&) {
        if (active.fetch_add(1) != 0) { ++overlaps; }
        std::this_thread::yield();
        ++delivered;
        --active;
    };
    {
        std::vector<std::jthread> workers;
        for (int i = 0; i < 4; ++i) {
            workers.emplace_back([&] {
                for (int message = 0; message < 20; ++message) { add_report(p, 3, "worker progress"); }
            });
        }
    }
    EXPECT_EQ(delivered, 80);
    EXPECT_EQ(overlaps, 0);
}

TEST(Feedback, LoadingIdentityRecordsIsAnExplicitlyReportedWrite)
{
    run_session session(":memory:");
    auto& p = session.parameters();
    p.active_operation = operation::load_inputs;
    std::vector<progress_event> events;
    p.progress = [&](const progress_event& event) { events.push_back(event); };
    execute(session.database(), R"(
        CREATE TABLE solve_settings (
            SOLVE_SETTING_ID INTEGER PRIMARY KEY, ALL_EXP_FITTED TEXT,
            PARAMETERS_SOLVED_FOR TEXT, REACTIONS_ENABLED TEXT,
            SCATTER_METHOD TEXT, X_RESOLUTION INTEGER, Z_RESOLUTION INTEGER);
        CREATE TABLE solutions (SOLUTION_ID INTEGER PRIMARY KEY,
            SOLVE_SETTING_ID INTEGER, EXPERIMENT_NAME TEXT, INLET_COND_ID INTEGER);
    )");
    p.experiment_runs.emplace_back();
    p.experiment_runs[0].name = "fixture";
    for (int inlet = 1; inlet <= 4; ++inlet) {
        p.experiments.emplace_back();
        p.experiments.back().run = &p.experiment_runs[0];
        p.experiments.back().INLET_COND_ID = inlet;
    }
    get_solve_settings_ID_from_db(p, session.database());
    get_SOLUTION_IDs_from_db(p, session.database());
    EXPECT_GT(p.SOLVE_SETTING_ID, 0);
    EXPECT_GT(p.experiments[0].SOLUTION_ID, 0);
    EXPECT_DOUBLE_EQ(scalar(session.database(), "SELECT count(*) FROM solve_settings;"), 1.0);
    EXPECT_DOUBLE_EQ(scalar(session.database(), "SELECT count(*) FROM solutions;"), 4.0);
    const auto creations = std::count_if(events.begin(), events.end(), [](const progress_event& event) {
        return event.action == operation::load_inputs && event.message.find("Creating ") == 0;
    });
    EXPECT_EQ(creations, 2); // One settings message and one summary for all four profiles.
    EXPECT_EQ(std::count_if(events.begin(), events.end(), [](const progress_event& event) {
        return event.message == "Creating 4 solution records";
    }), 1);
    events.clear();
    get_solve_settings_ID_from_db(p, session.database());
    get_SOLUTION_IDs_from_db(p, session.database());
    EXPECT_FALSE(std::any_of(events.begin(), events.end(), [](const progress_event& event) {
        return event.message.find("Creating ") == 0;
    }));
}

TEST(Feedback, DatabaseWriteFailuresHaveSqliteCodes)
{
    run_session session(":memory:");
    auto& p = session.parameters();
    p.active_operation = operation::save_model_profiles;
    try {
        delete_values_from_db(p, session.database(), "missing_table", "1=1");
        FAIL() << "A failed write must not be silently accepted";
    } catch (const workflow_error& error) {
        EXPECT_EQ(error.code, error_code::database);
        EXPECT_EQ(error.action, operation::save_model_profiles);
        EXPECT_EQ(error.sqlite_code, SQLITE_ERROR);
    }
}
} // namespace
