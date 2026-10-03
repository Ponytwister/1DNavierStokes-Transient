#include <background_runner.h>

namespace tsensor_workflow {

background_runner::~background_runner()
{
    request_cancel();
    wait();
}

void background_runner::start(const std::filesystem::path& database, progress_callback observer, std::optional<control_values> controls)
{
    // Resolve before dispatch so a later process working-directory change cannot
    // change which database the worker opens.
    if (database.empty()) {
        throw workflow_error(operation::open_database, error_code::invalid_input,
                             "Database path must not be empty");
    }
    const auto path = std::filesystem::absolute(database).lexically_normal();
    std::lock_guard lock(mutex_);
    if (state_ != background_state::idle) {
        throw workflow_error(operation::none, error_code::invalid_state,
                             "Take the previous outcome before starting another run");
    }
    cancellation_ = std::stop_source{};
    events_.clear();
    state_ = background_state::running;
    try {
        worker_ = std::jthread([this, path, observer = std::move(observer), controls = std::move(controls),
                               token = cancellation_.get_token()] {
            background_result outcome;
            auto state = background_state::completed;
            try {
                auto progress = [this, observer](const progress_event& event) mutable {
                    {
                        std::lock_guard lock(mutex_);
                        // Bounded recent history; final outcome is stored separately.
                        if (events_.size() == 512) events_.pop_front();
                        events_.push_back(event);
                    }
                    if (observer) observer(event);
                };
                auto session = std::make_unique<run_session>(path, std::move(progress), token);
                session->load_inputs(controls);
                outcome.result = session->run();
                // The session can outlive this runner; remove its capturing callback.
                session->parameters().progress = {};
                outcome.session = std::move(session);
            } catch (const workflow_error& error) {
                state = error.code == error_code::cancelled
                      ? background_state::cancelled : background_state::failed;
                outcome.failure = std::current_exception();
            } catch (...) {
                state = background_state::failed;
                outcome.failure = std::current_exception();
            }
            std::lock_guard lock(mutex_);
            // A cancellation accepted before completion wins even if it arrived
            // after the last solver checkpoint. Never expose a partial session.
            if (state == background_state::completed && token.stop_requested()) {
                state = background_state::cancelled;
                outcome = {};
                outcome.failure = std::make_exception_ptr(workflow_error(
                    operation::solve, error_code::cancelled, "Run cancelled"));
            }
            if (outcome.session) outcome.session->parameters().cancellation = {};
            outcome_ = std::move(outcome);
            state_ = state;
        });
    } catch (...) {
        state_ = background_state::idle;
        throw;
    }
}

bool background_runner::request_cancel()
{
    std::lock_guard lock(mutex_);
    return state_ == background_state::running && cancellation_.request_stop();
}

background_state background_runner::status() const
{
    std::lock_guard lock(mutex_);
    return state_;
}

std::vector<progress_event> background_runner::drain_events()
{
    std::lock_guard lock(mutex_);
    std::vector<progress_event> events(events_.begin(), events_.end());
    events_.clear();
    return events;
}

void background_runner::wait()
{
    if (worker_.joinable()) worker_.join();
}

background_result background_runner::take_result()
{
    {
        std::lock_guard lock(mutex_);
        if (state_ == background_state::running || state_ == background_state::idle) {
            throw workflow_error(operation::none, error_code::invalid_state,
                                 "No finished background outcome is available");
        }
    }
    wait();
    std::lock_guard lock(mutex_);
    auto outcome = std::move(outcome_);
    outcome_ = {};
    state_ = background_state::idle;
    return outcome;
}

} // namespace tsensor_workflow
