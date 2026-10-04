# Windows Qt kit selection

The initial desktop UI will use the following fixed kit, recorded in
[`cmake/QtKit.cmake`](../cmake/QtKit.cmake):

| Item | Selection |
| --- | --- |
| Qt | 6.10.3, official Windows MinGW 64-bit binary kit |
| Compiler | Qt-distributed MinGW-w64 GCC 13.1.0, x86_64 |
| UI module | Qt Widgets (with its Core and Gui dependencies) |
| Linking | Shared Qt libraries |
| Project language | GNU C++20, extensions enabled |
| Build system | Existing CMake 4.2+ and MinGW Makefiles |

Qt 6.10 documents MinGW-w64 13.1 support on Windows 10 (1809+) and Windows 11.
Qt Widgets supplies the desktop controls needed for file selection, buttons,
logs, and result tables. No Qt Quick, WebEngine, Charts, or Qt SQL dependency is
needed for the initial UI; the core keeps its existing SQLite library.

Sources:

- [Qt 6.10 Windows compiler support](https://doc.qt.io/qt-6.10/windows.html)
- [Qt 6.10.3 Widgets, CMake integration, and licensing](https://doc.qt.io/qt-6.10/qtwidgets-index.html)
- [Qt 6.10.3 release archive](https://download.qt.io/archive/qt/6.10/6.10.3/)

## Current machine and installation

The selected kit is installed side by side at the paths below. Executable checks
report Qt `6.10.3`, GCC `13.1.0`, and target `x86_64-w64-mingw32`.
The existing MSYS2 GCC 12.1.0 (Rev3) at `C:\msys64\mingw64\bin\g++.exe`
remains available for the original terminal presets. System-wide PATH was not changed.

Installation used the pinned `aqtinstall==3.3.0` helper in an isolated Python 3.12
virtual environment under ignored `out/qt-setup/venv`. It downloaded Qt's public
archives; the helper is not an application dependency. Equivalent commands after
creating a Python 3.10+ virtual environment and installing that helper are:

```powershell
python -m aqt install-qt --outputdir C:\Qt windows desktop 6.10.3 win64_mingw --archives qtbase d3dcompiler_47 opengl32sw
python -m aqt install-tool --outputdir C:\Qt windows desktop tools_mingw1310
```

The Qt Base archive includes Core, Gui, Widgets, their build tools, and plugins.
The graphics runtime archives are included for Qt's Windows graphics support.
[aqtinstall documentation](https://aqtinstall.readthedocs.io/en/v3.3.0/cli.html)
describes these installation commands. The official Qt installer remains an
alternative:

Install the **Qt 6.10.3 MinGW 64-bit desktop component** and the matching
**MinGW 13.1.0 64-bit compiler tools** through the Qt installer. Use a side-by-side
installation, for example:

```text
C:\Qt\6.10.3\mingw_64
C:\Qt\Tools\mingw1310_64
```

If that exact release is not offered in the default installer view, check its
archived-version selection. Do not silently substitute another Qt or compiler
version; revise the kit pin explicitly if a different kit is needed. Qt Creator
is optional; CMake can build the application directly.

Use the compiler supplied with this kit for the GUI executable **and all its
linked project libraries**, including ALGLIB and SQLite. Do not link the new UI
against archives from `out/codex-debug`, `out/codex-release`, or `out/build`.
Do not mix the old MSYS2 runtime DLLs into the new application's PATH or package.
No system-wide PATH change or MSYS2 upgrade is required.

After installation, check the actual paths and versions in PowerShell:

```powershell
& 'C:\Qt\Tools\mingw1310_64\bin\g++.exe' -dumpfullversion
& 'C:\Qt\Tools\mingw1310_64\bin\g++.exe' -dumpmachine
& 'C:\Qt\6.10.3\mingw_64\bin\qmake.exe' -query QT_VERSION
& 'C:\Qt\6.10.3\mingw_64\bin\qmake.exe' -query QT_INSTALL_PREFIX
```

Expected versions are `13.1.0` and `6.10.3`, with an `x86_64-w64-mingw32`
compiler target and the selected Qt installation prefix.

## Building and launching

Qt remains optional: `NAVIER_BUILD_GUI` defaults to `OFF`. The original
`mingw-debug` and `mingw-release` presets require no Qt. The GUI presets enable
it, enforce the pinned Qt version and GNU compiler/architecture, and rebuild all
project libraries in fresh `out/qt-debug` and `out/qt-release` directories.

From the repository root:

```powershell
cmake --preset qt-mingw-debug
cmake --build --preset qt-mingw-debug --parallel 2
ctest --preset qt-mingw-debug
```

Use `qt-mingw-release` for the release equivalents. For offline configuration,
append the existing `-DFETCHCONTENT_SOURCE_DIR_GOOGLETEST=...` override described
in the main README. On this machine it points at
`D:/GitHub/1DNavierStokes-Transient/out/build/_deps/googletest-src`.

The GUI executable is `out/qt-debug/gui/NavierGui.exe` (or the release equivalent).
The terminal executable remains `Navier.exe` at the build-directory root. CMake
presets set the kit PATH and Qt plugin directory for configure/build/test only.
For a manual development launch in PowerShell:

```powershell
$env:PATH = "C:\Qt\Tools\mingw1310_64\bin;C:\Qt\6.10.3\mingw_64\bin;" + $env:PATH
$env:QT_PLUGIN_PATH = "C:\Qt\6.10.3\mingw_64\plugins"
& .\out\qt-debug\gui\NavierGui.exe
```

## Using the desktop window

Choose an existing experiment database, then click **Run**. The window starts
with no database selected and never opens the repository database automatically.
The output directory can be selected or entered; export creates it if necessary.
Inputs are locked while a calculation or save is active. The progress indicator
shows activity, not a completion percentage; completed model evaluations are
reported separately. The log keeps the most recent 500 lines.

The result table lists parameter source, name, and value. The summary reports
optimizer iterations and the termination code; "finished" does not imply a
converged fit. The result database is displayed above the table. Choosing another
database clears the old results to prevent saving them to the wrong destination.
Starting a new run also replaces the previous result session.

Use **File > Open setup...** and **File > Save setup...** for setup files.
**File > Database** and **File > Output directory** provide editable paths and
Browse buttons. These actions are locked while work or control editing is active.

Before running, open the **Model controls** tab to edit the selected database's
`model_controls` rows. Boolean controls use true/false choices; resolution,
padding, iteration limits, report level and convergence tolerance are validated
when saving. The experiment dropdown lists names from the database experiments
table with checkboxes for multiple selections. At least one experiment is required;
saved names missing from the table must be deselected. Names containing whitespace
cannot be selected because the model uses a space-separated list. Existing selection
order is preserved, with new selections appended. The optional global-parameter dropdown offers `p1`,
`kon1`, `keq1`, `left_edge`, `width`, and `QE1`. If any experiment in the database
has more than one space-separated reaction in `REACTIONS`, it also offers `p2`,
`kon2`, `keq2`, and `QE2`, regardless of which experiments are selected. Clearing
all parameter selections stores NULL. Parameter-link semantics are unchanged. Each control has a dedicated form row, with checkboxes for
true/false values. Hover over a control for a description. Only
`universal_solve_for` may be blank (stored as SQL NULL); all other values are
required. Both apply actions validate every displayed field, including unchanged
values. `debug_level` accepts integers 0 through 6, and `convergence_epsx` must be
finite and strictly between 0 and 1e-3. These are editor constraints; model equations
and CLI loading behavior are unchanged. Validation does not establish whether a
combination is physically appropriate or will converge.

Use **Use values**, **Update default**, or **Cancel** to finish editing and return
to **Run and results**. Run and file actions remain disabled until editing finishes.

Controls start from database defaults and remain in memory for subsequent runs.
**Use values** applies edits in memory; **Update default** also writes them to the
database; **Cancel** discards tab edits. Applying controls clears previous results to prevent
mixing results with newly edited inputs. The editor is unavailable while running
or saving results. Reads and writes run in a worker while the tab stays
responsive. Close/Cancel waits until a pending database operation finishes.

Both the original `Parameter`/`Setting` layout and the test `criterion`/`value`
layout are supported, following the loader's first-two-column convention.
No schema migration is needed. The editor updates existing recognized rows only;
unknown rows, unchanged values, extra columns, and experiment/result tables are preserved.
Unknown rows are not displayed as controls. The legacy `save_normalized_profiles`
and `save_model_profiles` rows are preserved in existing databases but are not
editable controls or applied by desktop runs. **Save profiles** remains available
after a successful run and saves only when explicitly clicked, regardless of the
legacy flag. CLI behavior is unchanged.
Updates use one transaction and roll back on failure. If another application
changes the controls after loading, saving refuses to overwrite them; close and
reopen the editor to reload. No writes occur when merely opening or cancelling.

The loader now honors `run_solver=false`: it evaluates the model without fitting,
using the existing non-optimizer workflow. Previously false left the default true
value unchanged. Equations, numerical units and parameter-link semantics are
unchanged.

After a successful run:

- **Export report** writes the existing text report into the output directory,
  replacing a report with the same experiment-derived filename if present.
- **Save profiles** writes result profiles into the run's database. This button
  is disabled when the database's `save_model_profiles` setting is false.
- **Save fitted inputs...** asks for confirmation before replacing fitted initial
  inputs in that database. It is separate from saving result profiles.

Export and saves run off the UI thread and are serialized with calculations.
Errors appear in the status and log; a failed export/save keeps results available
for retry. Nothing is exported or saved automatically after a calculation.
Loading can still create identity records, as documented for the core workflow.
Existing partial-write behavior is unchanged: a failing database save does not
roll back statements already committed.

**Cancel** requests cooperative cancellation. Closing during a calculation
requests cancellation and defers closing until workers finish; closing during a
save lets that save finish. If saving fails, the window stays open with the error
and results available for retry; close again to exit without retrying.
The event loop continues in both cases. Cancellation
latency depends on the existing core checkpoints. Fitted-parameter input editing
and plotting are not yet included.

`Application.GuiStartup` tests window construction and the event loop with Qt's
offscreen plugin. GUI integration tests use disposable synthetic databases
and exercise controls, explicit persistence, cancel/close, invalid-input recovery,
export retries (including failure during closing), event-loop cancellation,
and the disabled-profile-save setting. These run alongside the
39 existing component/workflow tests. The same numerical coverage limitations
apply; these checks do not validate general parameter recovery or deployment.
Controls-editor tests additionally cover explicit save/cancel, invalid values,
the real column-name layout, NULL preservation, transaction rollback, conflicting
external edits, and running without fitting after editing the switch.

If installing outside `C:\Qt`, use local preset overrides for the C/C++ compiler,
make tool, Qt prefix, PATH, and QT_PLUGIN_PATH together. Keep versions consistent
with `cmake/QtKit.cmake`; do not change compilers inside an existing build tree.

## Licensing and redistribution

Qt Widgets offers commercial, LGPLv3, and GPLv2 licensing options. The planned
shared-library setup supports using the LGPLv3 distribution route; this does not
assign a license to this repository. Distribution must include the applicable Qt
licenses/notices and satisfy the chosen license's source and replacement/relinking
requirements. Record the actual distributed Qt components and their third-party
notices when packaging. Static linking is not part of this kit decision.

## Portable Windows package

From the repository root, build the Release ZIP with the pinned kit:

```powershell
cmake --preset qt-mingw-release
cmake --build --preset qt-mingw-release --parallel 2
ctest --preset qt-mingw-release
cmake --build --preset qt-mingw-release --target package
./tools/test-portable-package.ps1 -Archive ./out/qt-release/packages/NavierGui-0.1.0-windows-x64-Release.zip
```

The package target uses CMake's CPack directly, avoiding PowerShell's possible
Chocolatey `cpack` alias. Version `0.1.0` identifies this initial development
package format, not a stable application release. All packaging output stays in
the ignored Qt build directory. Reconfigure first after changing install rules.

Extract the entire ZIP and launch `bin/NavierGui.exe`. The package contains the
GUI, shared Qt runtime, Windows platform/style/image plugins, `qt.conf`, the
matching compiler runtimes, dependency records and a usage README. No database
or experiment output is included. The CLI remains available in the build tree.
Qt's deployment helper determines runtime files; compiler DLLs are copied from
the selected compiler directory rather than whichever compiler is on PATH.
Unused network/touch plugin groups and translations are excluded for this UI.

The smoke script extracts into a unique folder beside the ZIP, clears Qt/QML
environment overrides and limits PATH to Windows system folders. It launches
the actual Windows platform plugin from a different working directory, checks
successful exit within 15 seconds, and restores the environment. It retains
the extraction for inspection. It never selects or opens a database.
This is startup/relocation validation on the development machine, not a clean-VM
test, full packaged solver validation, or a signed installer.

Read [the packaged README](portable-package.md) before redistribution. In
particular, the existing ALGLIB sources declare GPL v2 or later, and this
repository has no application license selected. The package includes available
compiler/Eigen/ALGLIB notices and Qt's supplied SBOM; completing corresponding
source and Qt/third-party license materials remains a public-release task.

Deployment API reference:
[Qt 6.10 deployment script](https://doc.qt.io/qt-6.10/qt-generate-deploy-app-script.html).

### Experiments tab

The **Experiments** tab displays all rows and columns of the selected database's
`experiments` table, ordered by `NAME`. SQL NULL is shown explicitly. The table
refreshes when opened or when the database changes; a read error clears stale
rows. Browsing experiments does not load a model or create solution records.

Channel dimensions can be edited and saved per experiment after applying the
[channel dimensions migration](channel-dimensions.md). Dimensions are in meters
and may differ between experiments in a combined run.

### Reactions and Species tabs

The **Reactions** and **Species** tabs display all columns of the selected
database's `reactions` and `species` tables, ordered by `REACTION_NAME` and
`SPECIES_NAME`, respectively. Table cells are read-only; use **Add...** or select one row and choose **Modify...**
to open the editor. Values, column names, and stored units are not transformed. SQL NULL displays as `NULL` with a tooltip.

Opening a tab, selecting another database, or clicking **Refresh** reloads the
view. A missing database/table or read error clears old rows and shows the error;
an empty table retains its column headers. Reads open existing databases only and
close their connections afterward. They do not create solution records or change
reaction rates, species properties, parameter links, or unsaved control edits.
Both tabs and their Refresh buttons are disabled while a run or result save is
active, following the main window's existing work-state controls.

The reference-table editors save only when **Save** is pressed; **Cancel** discards
edits. Names are fixed when modifying a row so existing experiment/reaction
references remain intact. New names must be unique, nonempty, and contain no
whitespace or apostrophes (the model uses space-separated names and SQL lookups).
Species type must be `molecule` or `particle`; numeric fields require finite
numbers or explicit SQL NULL. Reaction species must exist in the database;
coefficients and exponents require one number per species, and `Ks` requires two
numbers, matching the loader's forward/reverse rate array. These checks do not
establish physical suitability or convergence.

Writes use bound parameters and a transaction. Modify compares the complete
original row, including unknown columns, to detect concurrent changes; on conflict,
cancel and Refresh before retrying. Extra columns are retained on Modify and use
database defaults on Add. A successful save clears previous run results. Add and
Modify are unavailable while model controls or channel dimensions have pending
edits, or while a calculation/result save is active. No schema migration is needed.
