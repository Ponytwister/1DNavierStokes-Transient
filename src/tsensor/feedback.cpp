#include <tsensor.h>

namespace tsensor_workflow {

const char* operation_name(operation action) noexcept
{
    switch (action) {
        case operation::open_database: return "Open database";
        case operation::load_inputs: return "Load inputs";
        case operation::solve: return "Run model/fit";
        case operation::export_results: return "Export results";
        case operation::save_model_profiles: return "Save model profiles";
        case operation::save_fitted_parameters: return "Save fitted parameters";
        default: return "Model operation";
    }
}

} // namespace tsensor_workflow

void publish_event(parameters_t& p, tsensor_workflow::progress_event event)
{
    std::lock_guard<std::recursive_mutex> lock(p.report_mutex);
    if (!p.progress) { return; }
    try {
        p.progress(event);
    } catch (...) {
        // A display failure must not change solver or persistence behavior.
        // The caller can inspect the failure after the synchronous operation.
        p.progress_failure = std::current_exception();
        p.progress = {};
    }
}

void check_cancellation(const parameters_t& p)
{
    if (p.cancellation.stop_requested()) {
        throw tsensor_workflow::workflow_error(p.active_operation,
            tsensor_workflow::error_code::cancelled, "Run cancelled");
    }
}
