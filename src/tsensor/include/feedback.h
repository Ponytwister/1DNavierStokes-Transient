#pragma once

#include <functional>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

namespace tsensor_workflow {

enum class operation { none, open_database, load_inputs, solve, export_results,
                       save_model_profiles, save_fitted_parameters };
enum class event_kind { started, message, evaluation, completed, failed };
enum class error_code { invalid_state, invalid_input, database, solver, io, internal };

struct progress_event {
    event_kind kind;
    operation action;
    std::string message;
    int detail_level = 3;
    // Completed residual evaluations, not optimizer iterations or a percentage.
    std::optional<int> evaluations;
    std::optional<error_code> error;
    std::optional<int> sqlite_code;
};

using progress_callback = std::function<void(const progress_event&)>;

class workflow_error : public std::runtime_error {
public:
    workflow_error(operation action, error_code code, std::string message,
                   std::optional<int> sqlite_code = {})
        : std::runtime_error(std::move(message)), action(action), code(code), sqlite_code(sqlite_code) {}
    const operation action;
    const error_code code;
    const std::optional<int> sqlite_code;
};

const char* operation_name(operation action) noexcept;

} // namespace tsensor_workflow
