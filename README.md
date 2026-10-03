# 1D Navier–Stokes transient model

C++ numerical modeling and parameter fitting using ALGLIB and SQLite. The application
reads experiment settings and profiles from `navier.db`, runs the configured model/fit,
and exports results.

## Build and test on Windows

Prerequisites: CMake 4.2 or later, Git (for the first GoogleTest download), and an MSYS2
MinGW64 toolchain providing `gcc`, `g++`, and `mingw32-make`. Put the toolchain's `bin`
directory on PATH, including when running tests so Windows can find its runtime DLLs.
For the usual installation, this PowerShell command updates only the current session:

```powershell
$env:PATH = "C:\msys64\mingw64\bin;" + $env:PATH
cmake --preset mingw-debug
cmake --build --preset mingw-debug --parallel 2
ctest --preset mingw-debug
```

Run these commands from the repository root. Use `mingw-release` in all three commands
to build and test the optimized configuration. Presets use separate `out/codex-debug`
and `out/codex-release` directories, leaving the existing `out/build` cache alone.
GNU C++20 extensions are enabled: the source currently uses GCC-specific headers,
variable-length arrays, and numeric literals. MSVC compatibility is not claimed.

The optional desktop UI uses Qt 6.10.3 Widgets with Qt's MinGW 13.1.0 64-bit kit.
See [Qt kit setup](docs/qt-setup.md) for installation, build/launch commands, and
licensing. Use the `qt-mingw-debug` or `qt-mingw-release` presets to build the
`NavierGui` application alongside the terminal executable. It provides database
and output selection, Run/Cancel, progress, parameter results, and explicit saves.
The **Model controls...** button edits controls held in memory, initially loaded
from the selected database. **Use values** applies edits for subsequent runs without
writing defaults; **Update default** also saves them to the database; **Cancel**
discards dialog edits. Changing the database resets the active controls.
The original presets remain terminal-only and require no Qt.

For a portable Windows ZIP, follow the [packaging instructions](docs/qt-setup.md#portable-windows-package).

GoogleTest is pinned to the existing `release-1.11.0` tag. Configuration normally fetches
it from GitHub. For offline configuration, supply an existing checkout of that version:

```powershell
cmake --preset mingw-debug -DFETCHCONTENT_SOURCE_DIR_GOOGLETEST="D:/path/to/googletest"
```

This machine already has a checkout under `out/build/_deps/googletest-src`; it can be
used as the override. The override is local CMake cache state, not part of the shared
preset. To build only the application without fetching GoogleTest, configure with
`-DBUILD_TESTING=OFF`; reconfigure with `-DBUILD_TESTING=ON` before running tests.

CTest discovers the individual GoogleTest cases. A failing assertion or no discovered
tests causes a failing test command. CLI checks cover help and invalid-input paths,
plus a full synthetic solve with both save-prompt answers in a disposable database.

## Test coverage and numerical fixtures

The initial suite checks linked parameter writes/unlinking, cycle rejection, interpolation
and its domain, molecular unit conversion, inlet initialization, and SQLite control loading.
See `tests/fixtures/README.md` for input formats, expected-value derivations, and tolerances.
The original cases are component regressions. The workflow fixture below adds a
complete transient solve in a deliberately limited uniform-equilibrium case.
Path tests cover option parsing, exports to disposable directories, file-open errors,
and CLI rejection of missing database files without creating an empty database.
Session tests cover independent state and links, owned concentration buffers,
partial-load failures, SQLite statement/connection cleanup, and worker exceptions.
Feedback tests cover silent core execution, lifecycle/error events, serialized
callbacks, evaluation counts, returned parameter snapshots, and explicit saves.

`workflow.sql` adds synchronous and background load/solve/export/save/reload checks,
cancellation after a real residual evaluation, and recovery from invalid inputs.
It checks the independent uniform-concentration solution, export contents, persisted
profiles, explicit fitted-input saving, replacement rather than duplication on
repeated saves, and preservation of unrelated experiment records. The CLI runs from
a different directory with explicit paths and answers both `n` and `y` to saving.
All test databases are synthetic and disposable; tests never run against `navier.db`.

This is not general solver validation: the case has no concentration gradient or
net reaction, and its fitted parameter is intentionally unidentifiable. It does not
establish accuracy for diffusion fronts, reacting systems, scattering corrections,
zero-profile comparisons, or recovery of parameters from experimental data. No
performance benchmark or sanitizer run is implied.

## Running the application

In a terminal, the model evaluation counter updates on one line. Redirected
output records only the final count, keeping log files free of repeated updates.

Pass `--database PATH` and `--output-dir DIRECTORY` to choose the existing experiment
database and export location. For example, from the repository root:

```powershell
.\out\codex-debug\Navier.exe --database ".\navier.db" --output-dir ".\out\results"
```

Relative paths are resolved against the launch working directory. Absolute paths
allow launching from any directory, including through a desktop shortcut. Quote
paths containing spaces. Use `--help` (or `-h`) to display usage without opening a
database. Unknown options and missing path values exit with an error.

With no options, the legacy defaults remain `../../navier.db` and `..`; launch from
`out/codex-debug` or `out/codex-release` to use them. Missing output directories are
created. Missing databases are rejected instead of silently creating empty files.
An output directory that cannot be created fails before model loading; export
open/write failures are also reported, though a failed write may leave a partial
file. Existing files with the same experiment-based name are still overwritten.
Experiment names containing `/`, `\`, or `:` are rejected at export so they cannot
redirect the file outside the selected directory.

Runs can write database results and prompt to save fitted parameters. Selecting a
different output directory does not isolate database writes: use a disposable copy
of the database when the original must be preserved. Automated checks do not run a
full solve or fit.

## Saving and opening a setup

In the desktop UI, choose a database and output directory, then select **Save
setup…** to write a `.navier.json` file containing both locations and a snapshot
of the active in-memory model controls. Saving a setup does not run the
model or change the database. Existing setup files are replaced atomically.

Use **Open setup…** to restore the paths and controls in memory without changing
database defaults. Incompatible control names and changed unknown controls are rejected.
Opening a setup clears previous run results; run again to calculate new results.
Setup actions are disabled while a run, results save, or controls editor is active.

The version 1 JSON format uses `format: "navier-setup"`, `version: 1`, `database`,
`outputDirectory`, and a `controls` array of `{ "name": "...", "value": "..." }`
rows. Values are strings or JSON `null` (SQL NULL). Saved paths are absolute;
manually supplied relative paths resolve against the setup file's directory.
The referenced database must still exist. A setup file does not contain experiment
data or results and is not a database backup. No database schema migration is
required. The terminal application's options are unchanged.

## Source layout

- `tsensor.cpp`: terminal entry point, database connection, and save prompt.
- `src/tsensor/application_options.cpp`: terminal path options and legacy defaults.
- `src/tsensor/workflow.cpp` and `include/workflow.h`: shared loading, fitting,
  export, save, and owned run sessions in the `tsensor` library.
- `src/tsensor/feedback.cpp` and `include/feedback.h`: progress delivery and typed
  workflow errors; terminal formatting lives in `tsensor.cpp`.
- `src/tsensor/include/tsensor.h`: parameters, experiment structures, and interfaces.
- `src/tsensor/sqlite_interface.cpp`: database I/O, numerical helpers, and model routines.
- `src/alglib-cpp`, `src/sqlite3`, `src/eigen-3.4.0`: bundled dependencies.
- `tests`: GoogleTest suite and small versioned fixtures.
- `in`, `navier.db`, `out`: experiment inputs, working database, and generated results.

Generated files are already tracked in parts of this repository. Stage source changes
selectively; repository artifact cleanup is separate from this setup.

## Shared application workflow

`tsensor_workflow::run_session` owns one database connection and one run's
`parameters_t`, optimizer buffers, and per-experiment concentration buffers.
A terminal or future UI caller can use the same operations:

```cpp
tsensor_workflow::run_session session(database_path);
session.load_inputs();
const auto result = session.run();
const auto written_file = session.export_results(output_directory);
session.save_model_profiles();
// Call session.save_fitted_parameters() only when the user chooses to save.
```

Destruction frees the buffers and closes the database on success or exception;
there is no manual cleanup call. `save_model_profiles` still respects the loaded
`save_model_profiles` setting. Loading and saving can write database records;
cleanup does not undo already committed writes or exported files.

Use a **fresh session for each new run**. A session accepts one load followed by
one solve, then exports/saves. A failed load/solve marks it failed; further
load/solve/save attempts are rejected. Export or save failures leave completed
results available for retry. Session and parameter-state objects cannot be copied
or moved because parameter links refer to objects inside them. Borrowed parameter
references/database handles must not outlive the session; do not resize linked
containers or change parameters while model workers are running.

The low-level component functions take an explicit `parameters_t&`; there is no
global parameter object. SQLite row callbacks and ALGLIB receive that same state
through their context pointers. SQLite callback exceptions are captured and
rethrown after SQLite finalizes the active statement; error-message buffers are
also owned. Existing model workers use joining thread owners, capture exceptions,
and report them back after joining the batch. Their shared report stack is
synchronized, and a zero hardware-concurrency hint falls back to one worker.

These preparation steps do not change numerical equations, units, boundary
conditions, optimizer settings, or convergence tolerances. Calls on one session
must be serialized; use the background runner below for UI work. Component
and lifetime tests do not establish full transient-solver or parameter-fit correctness.

## Progress, results, and errors

The core neither reads terminal input nor writes to stdout/stderr. With no
observer it is silent. Pass an optional callback to the session constructor:

```cpp
tsensor_workflow::run_session session(database_path,
    [](const tsensor_workflow::progress_event& event) {
        // Copy the event into a UI queue or record it in a log.
    });
```

Events identify their operation and kind: started, message, evaluation, completed,
or failed. Messages have the existing detail level; the model's debug setting
filters messages but does not suppress lifecycle/evaluation events. Evaluation
events count completed residual evaluations, **not** optimizer iterations or a
percentage of work. The total amount of fitting work is not known in advance.

Callbacks are synchronous and serialized per session. Some run on model workers;
a future GUI must marshal copies to its UI thread. Callbacks must be short and
must not reenter the session, mutate its state, or wait for its workers. If a
callback throws, it is disconnected and its exception is available through
`session.progress_failure()` after the operation returns. Display failures do not
change model results or database-write behavior. The terminal subscribes to these
events and retains the explicit save prompt.

`run()` returns an owned `run_result`: whether the optimizer ran, its iteration
count and ALGLIB termination code when applicable, completed residual evaluations,
and named parameter values from the returned optimizer vector. If optimization
is disabled, the existing branch is preserved and the parameter snapshot contains
its input vector. The summary can outlive the session; larger profile arrays remain
in the session's parameters. A completed operation is not a claim of convergence:
callers must inspect the termination code.

Session failures throw `workflow_error`, with `action`, `code`, and an optional
native `sqlite_code`, alongside the human-readable `what()` message. Categories
include invalid state/input, database, solver, I/O, internal failures, and cancellation. Failure
events carry the same fields. Callers should use these fields instead of parsing
message strings. SQL errors now follow one error path, including failed deletes
that previously only printed a warning and let saving continue.

## Database writes in the workflow

- Loading may insert missing `solve_settings` and `solutions` identity records.
  Progress messages identify these insertions. It does not save fitted inputs.
- `save_model_profiles()` respects its setting and writes `parameter_solutions`,
  updates `solutions` metrics, and replaces the selected `model_profile` rows.
- Only an explicit `save_fitted_parameters()` call updates fitted initial inputs
  in `alglib_input`, `species`, `reactions`, `experiments`, and `raw_profile`.

These are existing persistence semantics, now documented for future UI callers.
No transaction or rollback behavior is added: a failed operation may already have
committed earlier statements. Use a disposable database for automated checks.

## Background execution and cancellation

`background_runner` (in `<background_runner.h>`) owns one active load/solve task.
The terminal continues to use the synchronous API. A future GUI can start work,
poll `status()` and `drain_events()` on a timer, and call `request_cancel()` without
blocking its event loop. Do not call blocking `wait()` from that event loop.

```cpp
tsensor_workflow::background_runner runner;
runner.start(database_path);
// Later, from the UI timer:
const auto events = runner.drain_events(); // Owned copies; render on the UI thread.
const auto state = runner.status();
if (state != tsensor_workflow::background_state::running &&
    state != tsensor_workflow::background_state::idle) {
    auto outcome = runner.take_result(); // Joins the finished worker, resets to idle.
    if (outcome.failure) std::rethrow_exception(outcome.failure);
    // outcome.result is the summary; outcome.session owns the completed profiles.
    // Keep the session for explicit export/save actions after completion.
}
```

The worker resolves and copies the database path at start and exclusively owns
its session. No mutable parameters or database handles are exposed while it runs.
Disable editing controls in the eventual GUI while status is `running`; the runner
rejects another start until the previous outcome is taken. This does not lock out
other applications editing the same database. Starting, waiting, taking results,
and destruction belong to one controlling thread. Status, cancellation, and event
draining are safe from other threads. Keep one runner per UI workflow.

Events retain the most recent 512 entries, dropping older entries if the UI falls
behind. Terminal status and outcome remain available independently of this queue.
The optional observer has the same callback contract as the synchronous session;
if it throws, event delivery is disconnected and a successful returned session
retains `progress_failure()`. The completed session is detached from the runner's
callback and stop token, so it can safely outlive the runner.

Cancellation is cooperative. Checkpoints run before workflow operations and SQL
statements/row callbacks, between model batches, at time steps, and within diffusion
iterations. Workers join before cancellation is reported; ALGLIB's bundled callback
exception guard cleans up while propagating the cancellation exception. No equations,
optimizer settings, or floating-point convergence criteria change. A long SQL
statement, normalization pass, optimizer internal step, or callback can delay a
checkpoint; there is no fixed cancellation-latency guarantee. Destruction requests
cancellation and joins, so arrange UI shutdown accordingly.

A cancelled outcome has `background_state::cancelled` and a `workflow_error` with
`error_code::cancelled`, with no usable partial session or result. A request accepted
before successful completion is published wins that race. A simultaneous unrelated
failure remains a failure. The runner never exports or saves profiles/fitted inputs
automatically, on either success or cancellation. Loading may already have created
identity records; cancellation does not roll those writes back. A completed session
supports the existing explicit export/save operations. Synchronous callers may also
pass a `std::stop_token` as the third `run_session` constructor argument.

Tests cover cancellation, active-run exclusion, outcome retrieval/reuse, failure
propagation, destruction, and ALGLIB callback unwinding with an independent scalar
least-squares case. Step 6 also validates a full synthetic background workflow and
checks that a completed session can outlive its runner and still export/save.
See the numerical limitations under Test coverage above.
