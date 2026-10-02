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

The window is a minimal application shell. Database selection, Run/Cancel,
progress, and result controls are the next UI step. Starting it does not open a
database. `Application.GuiStartup` constructs the window, enters the event loop,
and exits automatically using Qt's offscreen plugin. The GUI kit also runs all
39 existing component/workflow tests. This smoke test is not visual inspection
or a standalone deployment check.

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

Qt runtime DLLs and the Windows platform plugin must accompany a distributed UI;
deployment validation is a later step. Qt source or binary packages will not be
vendored into the repository.
