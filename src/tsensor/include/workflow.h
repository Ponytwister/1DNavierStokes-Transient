#pragma once

#include <tsensor.h>

// Shared application operations, extracted from the terminal entry point.
// The caller supplies the existing global `p` and owns the open database.
// These operations retain console reporting and use global state: they are
// single-run, non-reentrant operations, not a thread-safe session API.
namespace tsensor_workflow {

// Load in the established order, including profile imports and solution-ID
// creation (which can write to the database), then normalize the profiles.
void load_inputs(sqlite3* db);

// Requires successfully loaded inputs. Retains the existing run_solver behavior
// and optimizer settings; ALGLIB exceptions propagate to the caller.
void run();

// Export the current results using the existing relative path/name convention.
void export_results();

// Save profiles only when p.save_model_profiles is enabled.
void save_model_profiles(sqlite3* db);

// Explicitly persist fitted parameters; never prompts for user input.
void save_fitted_parameters(sqlite3* db);

// Existing successful-run cleanup. Call once after the last export/save.
// This does not reset all global state or establish repeated-run support.
void release_run_resources();

} // namespace tsensor_workflow
