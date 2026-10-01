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
  Document changes to equations, units, boundary conditions, convergence criteria, or
  floating-point behavior, and validate them against analytical or independent references.
- For numerical changes, add a small independently justified regression case. Existing
  fixtures use `abs(actual - expected) <= 1e-12 + 1e-10 * abs(expected)` and require finite
  outputs. Explain any tolerance change instead of accepting newly generated results blindly.
- The initial suite covers helpers, inlet initialization, and control loading, not a full
  transient solve. Report that limitation when assessing solver-wide correctness.

## Change and validation policies

- Keep changes focused on the requested task. Avoid unrelated cleanup, renaming,
  formatting, and dependency upgrades; include supporting changes only when necessary.
- For bug fixes, add a regression test that demonstrates the failure before the fix when
  practical. Never weaken assertions or change expected results solely to make tests pass.
- Make database schema changes explicit and reproducible through migration scripts or
  equivalent versioned steps. Test migrations on disposable databases, verify preservation
  of existing experiment data, and document compatibility and recovery steps.
- When adding or upgrading a dependency, explain the need, pin its version, and document
  build and licensing implications. Avoid incidental changes to bundled dependencies.
- Report the checks actually run, their results, and remaining gaps. Distinguish component
  tests, full solver validation, and benchmarks; do not imply an unperformed check passed.
- Support performance claims with repeatable benchmarks using the same inputs, compiler
  configuration, and hardware. Report the method and measured results, and verify that
  the optimization preserves numerical correctness.

## Atomic commits

- When committing, make each commit one coherent, reviewable change with a single purpose.
  Keep the implementation, its relevant tests, and required documentation together.
- Split independent fixes, features, refactors, and formatting changes into separate commits.
  Do not split tightly coupled changes if doing so leaves an intermediate commit broken.
  Atomicity is about logical scope, not a fixed file count or line limit.
- Inspect the working tree and staged diff before every commit. Stage explicit paths or
  hunks so unrelated user edits and generated files are not swept into the commit.
  Commit pre-existing user work separately when the user requests it.
- Each code commit should build and pass the checks relevant to its scope. Validation
  must cover the state being committed, not rely on unrelated uncommitted changes.
  For documentation-only changes, review the diff and run `git diff --cached --check`.
- Write a commit message that states the change and its purpose. Report the commit hash
  and validation performed, including any checks that could not be completed.
- Treat local commits, pushing, merging, and history rewriting as distinct actions within
  the user's authorized scope. Permission to commit does not by itself authorize pushing
  or merging. These policies do not themselves authorize any of those actions.
  Do not amend, squash, or reorder existing commits unless requested.
