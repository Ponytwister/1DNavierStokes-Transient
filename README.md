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
tests causes a failing test command. CLI checks launch only help and invalid-input
paths; they do not run the interactive solver.

## Test coverage and numerical fixtures

The initial suite checks linked parameter writes/unlinking, cycle rejection, interpolation
and its domain, molecular unit conversion, inlet initialization, and SQLite control loading.
See `tests/fixtures/README.md` for input formats, expected-value derivations, and tolerances.
These are component regressions, not validation of a full transient solve or parameter fit.
Path tests cover option parsing, exports to disposable directories, file-open errors,
and CLI rejection of missing database files without creating an empty database.
Session tests cover independent state and links, owned concentration buffers,
partial-load failures, SQLite statement/connection cleanup, and worker exceptions.

## Running the application

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

## Source layout

- `tsensor.cpp`: terminal entry point, database connection, and save prompt.
- `src/tsensor/application_options.cpp`: terminal path options and legacy defaults.
- `src/tsensor/workflow.cpp` and `include/workflow.h`: shared loading, fitting,
  export, save, and owned run sessions in the `tsensor` library.
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
session.run();
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

This step changes ownership and error cleanup, not numerical equations, units,
boundary conditions, optimizer settings, or convergence tolerances. Console
reporting remains. Calls on one session must be serialized; background UI work,
cancellation, and structured progress reporting are later steps. Component and
lifetime tests do not establish full transient-solver or parameter-fit correctness.
