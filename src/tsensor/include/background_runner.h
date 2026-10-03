#pragma once

#include <workflow.h>
#include <thread>

namespace tsensor_workflow {

enum class background_state { idle, running, completed, cancelled, failed };

struct background_result {
    std::unique_ptr<run_session> session;
    std::optional<run_result> result;
    std::exception_ptr failure;
};

// start/wait/take_result/destruction belong to one controlling thread. status,
// request_cancel and drain_events may also be called from other threads.
// The worker exclusively owns the session until take_result joins it.
class background_runner {
public:
    ~background_runner();
    // Rejects an active run or an outcome not yet taken. The optional observer
    // has the same non-reentrant, short-running contract as run_session.
    void start(const std::filesystem::path& database, progress_callback observer = {},
               std::optional<control_values> controls = std::nullopt);
    bool request_cancel();
    background_state status() const;
    std::vector<progress_event> drain_events();
    void wait(); // Blocking; do not use on a GUI event loop while running.
    background_result take_result(); // Rejects running/idle; resets to idle.

private:
    mutable std::mutex mutex_;
    background_state state_ = background_state::idle;
    std::stop_source cancellation_;
    std::deque<progress_event> events_;
    background_result outcome_;
    std::jthread worker_; // Joined before any members it uses are destroyed.
};

} // namespace tsensor_workflow
