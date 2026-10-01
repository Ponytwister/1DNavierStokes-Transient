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
tests causes a failing test command. Tests do not launch the interactive application.

## Test coverage and numerical fixtures

The initial suite checks linked parameter writes/unlinking, cycle rejection, interpolation
and its domain, molecular unit conversion, inlet initialization, and SQLite control loading.
See `tests/fixtures/README.md` for input formats, expected-value derivations, and tolerances.
These are component regressions, not validation of a full transient solve or parameter fit.

## Running the application

`Navier.exe` currently opens `../../navier.db` relative to the working directory. With a
preset build, launch it from `out/codex-debug` or `out/codex-release`, not the repository
root. Runs can import profiles, write database results, produce files under `out`, and
prompt to save fitted parameters. Use a disposable copy of the project/data for experiments
when the original database must be preserved. There is no isolated command-line smoke
mode yet; automated verification uses the test executable instead.

## Source layout

- `tsensor.cpp`: terminal entry point, database connection, and save prompt.
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
`load_inputs`, `run`, `export_results`, and `save_model_profiles` in that order.
`save_fitted_parameters` is a separate, explicit action; the terminal application
still asks whether to perform it. After the last export/save, the existing
`release_run_resources` operation runs, and the caller closes the database.

This is an extraction of the current single-run workflow, not yet a GUI/session
API: global state, console reporting, relative paths, and existing error handling
remain. Loading can import profiles and create database records. The profile-save
operation respects `p.save_model_profiles`. Cleanup must only run once after a
successful load/run; it does not reset all state for another run. Callers must not
run these operations concurrently. Optimizer settings, numerical equations, and
the existing `run_solver` branches are unchanged. The unused callback context in
the direct-call branch is now explicitly null instead of indeterminate.
