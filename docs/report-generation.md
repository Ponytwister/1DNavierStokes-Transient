# Report generation from completed results

Choose **Prepare report...** after a run, or **Use current run** in Report.
No output directory is required to prepare the report. Generation runs in the
background and reads the completed session into an in-memory snapshot. It does
not rerun the solver, save database profiles, or create an intermediate report
file. Choose blocks, profile types, columns, and formatting, then use **Export
selected report...** or **Copy for Excel**. These use the same previewed data.

Selections for matching profile types, columns and blocks survive regeneration.
New types are selected by default. Formatting choices remain in memory. Starting
a new run or invalidating results clears a generated preview but keeps the
choices. Imported files remain independently usable. Settings are not persisted
across application restarts.

## Additional profiles and units

The new generation path separates two previously conflated rows:

- `Analytical_Zero_(<dye model units>)` exports `analytical_zero` unchanged on
  `channel_position`, with `window_size` samples. The solver stores it as the
  dimensionless analytical expression multiplied by `run.dye_conc`, which is in
  the dye's declared model units. The existing solver fills it only inside the
  central 20–80% channel interval and when `iterations > 1`; this report does not
  recompute it or fill the remaining initialized zeros.
- `Bound_Beads_(<bead input units>)` restores the previously commented formula:
  `species_out[Bound_Dye_1][x] / abs(reaction_1.coef[FITC]) * bead_unit`, where
  `bead_unit` is the existing conversion from bead model units to bead input units.
  It uses all `X` model cells. A zero or nonfinite coefficient is rejected.
  This is the original first-reaction/first-bound-species interpretation, not a
  new sum over multiple binding reactions.

Every `Species:<name>_(<model units>)` row exports `species_out[species][x]`
directly on the model grid, without a unit conversion. These include additional
bound species when present. The model-grid position is
`x * run.W * 1e6 / X` in micrometers, consistent with the existing report header.
Unbound and total bead labels in generated reports use their actual declared
input units instead of assuming wt%. No model equations or boundary conditions
change.

Generation preserves double round-trip precision (`max_digits10`, classic
locale) before optional UI formatting. Fixed/scientific decimal settings round
the exported presentation only. File import preserves the precision already
present in that file; it cannot restore omitted profiles or lost digits.

## Compatibility and validation

**Export legacy report**, the CLI, and `save_excel_output()` retain the old
13-row format and default numeric precision. In those files, the historical
`Bound_Beads_(wt%)` label still contains analytical zero for compatibility.
Importing such a file does not reinterpret or invent a true bound-bead profile.

An independent regression uses a bound-dye value of 12, a stoichiometric magnitude
of 3, bead model units mg/ml, and solution density 2 g/ml: 4 mg/ml bound beads
equals 0.2 wt%. It separately checks stored analytical-zero samples, full-grid
species values, total beads, no source mutation, no database writes, and legacy
output compatibility. GUI tests verify that generation requires a completed run,
needs no output folder, retains selections, and excludes deselected data from
the first exported file. These are component and synthetic-workflow checks, not
general solver validation.
