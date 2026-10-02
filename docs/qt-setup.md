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

Inspection found MSYS2 GCC 12.1.0 (Rev3) at
`C:\msys64\mingw64\bin\g++.exe`. No Qt installation was found in the checked
default locations or on PATH. That compiler remains the existing CLI baseline;
the new kit is selected from Qt's supported configurations, not yet locally
build-tested. The package-manager query was blocked by an MSYS signal-pipe error,
so the inspection also checked the package metadata and CMake directories directly.

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

## Integration boundary

This step records the dependency/toolchain selection. It does not install Qt,
add a GUI executable, or change the existing build presets. The kit file is not
yet enforced by the terminal build.

The next step adds optional GUI integration, using
`find_package(Qt6 ${NAVIER_QT_VERSION} EXACT REQUIRED COMPONENTS Widgets)` inside
the GUI option, and checks the GNU compiler version and 64-bit target. Separate
GUI presets must specify both C and C++ compiler paths, the matching make tool,
Qt's prefix, and fresh build directories. Rebuild the core and run the existing
39-test suite under that kit before treating local compatibility as verified.
Keep the terminal-only build usable without Qt installed.

## Licensing and redistribution

Qt Widgets offers commercial, LGPLv3, and GPLv2 licensing options. The planned
shared-library setup supports using the LGPLv3 distribution route; this does not
assign a license to this repository. Distribution must include the applicable Qt
licenses/notices and satisfy the chosen license's source and replacement/relinking
requirements. Record the actual distributed Qt components and their third-party
notices when packaging. Static linking is not part of this kit decision.

Qt runtime DLLs and the Windows platform plugin must accompany a distributed UI;
deployment validation is a later step. Qt source or binary packages will not be
vendored into the repository by this setup step.
