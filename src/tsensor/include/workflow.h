#pragma once

#include <tsensor.h>
#include <memory>

// Low-level operations take explicit state. For an application run, prefer
// run_session below, which owns state/connection and enforces operation order.
namespace tsensor_workflow {

struct parameter_value {
    std::string source;
    std::string name;
    double value;
};

struct run_result {
    bool optimizer_ran = false;
    int residual_evaluations = 0;
    std::optional<alglib::ae_int_t> optimizer_iterations;
    std::optional<alglib::ae_int_t> termination_type;
    // Owned snapshot of the returned optimizer vector (or inputs if disabled).
    std::vector<parameter_value> parameters;
};

// Load in the established order, including solve-settings/solution-ID creation
// (which can write to the database), then normalize the profiles.
void load_inputs(parameters_t& p, sqlite3* db);

// Requires successfully loaded inputs. Retains the existing run_solver behavior
// and optimizer settings; ALGLIB exceptions propagate to the caller.
run_result run(parameters_t& p);

// Export into an explicit directory (created if absent), retaining the existing
// experiment-based filename. Returns the written path; reports I/O errors by
// exception. Relative paths are relative to the caller's working directory.
std::filesystem::path export_results(parameters_t& p, const std::filesystem::path& output_directory);

// Save profiles only when p.save_model_profiles is enabled.
void save_model_profiles(parameters_t& p, sqlite3* db);

// Explicitly persist fitted parameters; never prompts for user input.
void save_fitted_parameters(parameters_t& p, sqlite3* db);

enum class session_state { empty, loaded, completed, failed };

// One load/solve per session. Construct a fresh session for another run, retaining
// the old one if its results are still needed. No copy/move: parameter links and
// model worker references must keep the same owner and addresses.
// Calls on the same session must be serialized; this is not a background UI API.
class run_session {
public:
    explicit run_session(const std::filesystem::path& database_path, progress_callback progress = {},
                         std::stop_token cancellation = {});
    ~run_session() = default;
    run_session(const run_session&) = delete;
    run_session& operator=(const run_session&) = delete;
    run_session(run_session&&) = delete;
    run_session& operator=(run_session&&) = delete;

    void load_inputs();
    run_result run();
    std::filesystem::path export_results(const std::filesystem::path& directory);
    void save_model_profiles();
    void save_fitted_parameters();

    session_state state() const noexcept { return state_; }
    // Inspect after an operation returns. A throwing observer is disconnected;
    // its exception is retained while the model operation continues normally.
    std::exception_ptr progress_failure() const noexcept { return parameters_.progress_failure; }
    parameters_t& parameters() noexcept { return parameters_; }
    const parameters_t& parameters() const noexcept { return parameters_; }
    // Borrowed for component operations. Do not close it; finalize any statements
    // before destroying the session. Database writes are not rolled back merely
    // because a later workflow operation fails.
    sqlite3* database() const noexcept { return database_.get(); }

private:
    struct close_database {
        void operator()(sqlite3* db) const noexcept { sqlite3_close_v2(db); }
    };
    void require_state(session_state expected, operation action) const;
    std::unique_ptr<sqlite3, close_database> database_;
    parameters_t parameters_{};
    session_state state_ = session_state::empty;
};

} // namespace tsensor_workflow
