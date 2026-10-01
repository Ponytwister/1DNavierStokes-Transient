# Working in this repository

- Build and test from the repository root using the MinGW presets:
  `cmake --preset mingw-debug`, `cmake --build --preset mingw-debug --parallel 2`,
  then `ctest --preset mingw-debug`. See README.md for prerequisites and offline setup.
- Use the release equivalents when validating optimization-sensitive numerical changes.
  The project currently requires GNU extensions; do not assume MSVC compatibility.
- Application entry point: `tsensor.cpp`. Model, numerical helpers, and database code:
  `src/tsensor/sqlite_interface.cpp`; structures and interfaces: `src/tsensor/include/tsensor.h`.
  Tests and small independent fixtures live in `tests/`.
- `src/alglib-cpp`, `src/sqlite3`, and `src/eigen-3.4.0` are bundled dependencies.
  Avoid editing or reformatting them unless the task specifically requires it.
- Preserve existing user edits. Build in the preset directories rather than reusing the
  tracked `out/build` tree. Do not stage generated binaries, caches, or experiment outputs.
- Tests must use in-memory or disposable databases. Do not run `Navier.exe` against
  `navier.db` as an automated check: application runs can change data and write outputs.
- Preserve numerical units, boundary conditions, and parameter-link semantics unless
  the requested change explicitly alters them. Device dimensions are in meters; other
  units are identified in model inputs. Do not infer an undocumented unit from a name.
- For numerical changes, add a small independently justified regression case. Existing
  fixtures use `abs(actual - expected) <= 1e-12 + 1e-10 * abs(expected)` and require finite
  outputs. Explain any tolerance change instead of accepting newly generated results blindly.
- The initial suite covers helpers, inlet initialization, and control loading, not a full
  transient solve. Report that limitation when assessing solver-wide correctness.
