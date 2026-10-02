# Navier desktop test package

Extract the entire ZIP to a writable folder and open `bin/NavierGui.exe`.
Keep its DLLs, `qt.conf`, and the plugins folder together. Qt and a compiler do
not need to be installed on the destination machine. This build targets x64
Windows 10 (1809+) / Windows 11 with Qt 6.10.3 and MinGW GCC 13.1.0.

Choose an existing experiment database and an output directory, then Run.
Use **Model controls...** before running to edit settings; **Save to database**
persists changes, while **Cancel** discards them. Saving controls clears old
results from the window so the next run uses the new settings.
No experiment database is shipped or selected automatically. Test with a copy
of your database: loading can create solution records. Report export, profile
saving, and fitted-input saving are explicit actions. Cancellation is cooperative;
closing waits for workers, and failed saves keep the window open for retry.

## Dependency and redistribution record

This is a local development/test package, not a completed public release.
The repository has not selected an application distribution license. Before
redistribution, resolve the combined application license, provide corresponding
source/build instructions and all applicable notices, and validate on a clean
Windows machine. No code signing or installer is provided.

- Qt 6.10.3 Core, Gui, Widgets and plugins are shared libraries. The included
  `notices/qt-sbom` records the installed Qt Base kit, including third-party
  components; it is a superset of the deployed binaries, not a package manifest.
  Obtain the matching Qt source, license texts and third-party notices from
  https://download.qt.io/archive/qt/6.10/6.10.3/ before redistribution. The SBOM
  is not a substitute for those texts or source obligations.
- ALGLIB 4.02.0 is statically linked. Its source headers declare GPL v2 or later;
  copies of GPL v2 and v3 are in `notices/alglib`. A Qt LGPL deployment alone
  does not settle the combined application's distribution requirements.
- Eigen 3.4.0 is header-only; its supplied license files are in `notices/eigen`.
- SQLite is statically linked; its source declares public-domain dedication.
  See https://sqlite.org/copyright.html.
- GCC 13.1.0 / MinGW-w64 runtime libraries are dynamically linked. Their supplied
  GCC runtime exception, GPL and MinGW/winpthreads notices are in
  `notices/compiler`.

The package intentionally contains no production databases, generated model
results, tests, compiler executables, or Qt development headers.
