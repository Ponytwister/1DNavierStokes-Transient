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
  export, save, and successful-run cleanup operations in the `tsensor` library.
- `src/tsensor/include/tsensor.h`: parameters, experiment structures, and interfaces.
- `src/tsensor/sqlite_interface.cpp`: database I/O, numerical helpers, and model routines.
- `src/alglib-cpp`, `src/sqlite3`, `src/eigen-3.4.0`: bundled dependencies.
- `tests`: GoogleTest suite and small versioned fixtures.
- `in`, `navier.db`, `out`: experiment inputs, working database, and generated results.

Generated files are already tracked in parts of this repository. Stage source changes
selectively; repository artifact cleanup is separate from this setup.

## Shared application workflow

`tsensor_workflow` separates the application operations from `main()`. A caller
defines the existing global `parameters_t p`, opens the database, then calls
`load_inputs`, `run`, `export_results(output_directory)`, and `save_model_profiles`
in that order. `export_results` accepts a filesystem path and returns the written
file path. A future UI can supply its selected directory directly; CLI parsing is
independent of the shared workflow.
`save_fitted_parameters` is a separate, explicit action; the terminal application
still asks whether to perform it. After the last export/save, the existing
`release_run_resources` operation runs, and the caller closes the database.

This is an extraction of the current single-run workflow, not yet a GUI/session
API: global state and console reporting remain, and general resource ownership is
still unchanged. The terminal reports path and standard exceptions, closes its
database on failure, and returns a nonzero exit code. Loading can create database
records. The profile-save
operation respects `p.save_model_profiles`. Cleanup must only run once after a
successful load/run; it does not reset all state for another run. Callers must not
run these operations concurrently. Optimizer settings, numerical equations, and
the existing `run_solver` branches are unchanged. The unused callback context in
the direct-call branch is now explicitly null instead of indeterminate.
