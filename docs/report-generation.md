# Report generation from completed results

Choose **Prepare report...** after a run, or **Use current run** in Report.
No output directory is required to prepare the report. Generation runs in the
background and reads the completed session into an in-memory snapshot. It does
not rerun the solver, save database profiles, or create an intermediate report
file. Choose blocks, profile types, columns, and formatting, then use **Export
selected report...** or **Copy for Excel**. These use the same previewed data.

Selections for matching profile types, columns and blocks survive regeneration.
New types are selected by default. Use **Legacy format** to select the historical
13-profile set and ordering, shared headers, and legacy labels. On the current
generated report, this preset produces the legacy presentation using expanded
report data. Both report generators now retain the same round-trip precision;
Fixed and Scientific remain optional presentation choices. **Save format...**
and **Load format...** store and
restore report selections in `.navier-report.json` files. Saving a setup can
associate one of these files with it; opening that setup loads the report format
for the next generated report. Formatting choices remain in memory. Starting a
new run or invalidating results clears a generated preview but keeps the
choices. Imported files remain independently usable. Save a format file or
associate one with a setup to persist choices across application restarts.

## Metadata columns and profile mapping

The former **D-A / integrals** column is named **Stacked metadata**. It remains
available alongside the six independently selectable metadata fields:

| Column | Value carried by this source row |
| --- | --- |
| `exp_DA` | Experimental_Derivative |
| `model_DA` | Numeric_Derivative |
| `exp_integral` | Experimental_Profile |
| `model_integral` | Numeric_Model_Profile |
| `analytic_exp_integral` | Experimental_Difference |
| `analytic_model_integral` | Numeric_Difference |

Each profile type has a dropdown that can assign one metadata field to that
profile row. Unselected profile types are greyed out and their dropdowns are
disabled. Each profile dropdown offers every metadata-column item except
**Stacked metadata**, including unchecked items. The dropdown and checkbox are
independent: an assignment places its value in **Stacked metadata**, while the
checkbox controls whether that field also appears as a separate column. When
that separate column is enabled, it shows the source field's values on every
selected profile row, regardless of the dropdown assignment. The dropdown
controls only which value appears in **Stacked metadata** for each profile.
Values are collected before profile filtering, so hiding a carrier row does not
remove metadata. Missing values remain `-`, never zero or values from another
block.

The **Stacked metadata** column remains independently selectable and
reorderable. The Legacy format keeps this historical stacked field enabled,
turns off the six individual metadata column checkboxes, and retains the default
profile-to-metadata assignments, preserving the legacy layout. Its units row
lists the assigned metadata names in the current metadata-column order,
including fields whose separate-column checkboxes are off. Profile labels with
run-name prefixes are matched to their metadata carrier type before values are
collected and placed.

Drag the items in **Metadata columns**, or use **Move up / Move down**, to set
their output order, including the **Stacked metadata** column. Saved report
formats retain this order and the profile mappings. The preview, clipboard,
export, regeneration, and imported reports use the same settings. Mapping only
rearranges existing values; it does not recalculate integrals or change units.

## Additional profiles and units

The new generation path separates two previously conflated rows:

- The bound-dye profile is labeled with the configured bound-dye species name,
  matching the free-dye profile's use of its configured species name.
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

Both expanded and legacy generators write doubles with `max_digits10` (17
significant digits for `double`, enough to round-trip the value). The expanded
generator uses the classic locale. **Original text** preserves source digits;
Fixed and Scientific options round the presentation to the chosen number of
decimal places. Formatting affects exported presentation only. File import
preserves the precision already present in that file; it cannot restore omitted
profiles or lost digits.

Expanded reports share one coordinate header for each distinct channel width,
like legacy reports. Each experiment's core profile set also shares one `sec`
axis/units row: model-grid profile values align with the shared coordinate
header, and sampled profiles align with the `sec` axis. Raw and rescaled source
profiles retain separate axis rows because their coordinate systems differ.

## Compatibility and validation

**Export legacy report**, the CLI, and `save_excel_output()` retain the old
13-row format and default numeric precision. Their analytical-zero row is
labeled `Analytical_Zero_(<dye model units>)`, matching the expanded Report tab.
The Report tab treats imported profile labels literally. In particular, the
Legacy format preset does not treat an older `Bound_Beads_(wt%)` row as the
analytical-zero profile; this removes the old report-label compatibility rule.
For the current generated report, the preset presents the input-unit total-bead
label as the legacy `Total_Beads_(wt%)` label. That alias is part of the saved
format settings and does not change the source profile or its values.

An independent regression uses a bound-dye value of 12, a stoichiometric magnitude
of 3, bead model units mg/ml, and solution density 2 g/ml: 4 mg/ml bound beads
equals 0.2 wt%. It separately checks stored analytical-zero samples, full-grid
species values, total beads, no source mutation, no database writes, and legacy
output compatibility. GUI tests verify that generation requires a completed run,
needs no output folder, retains selections, and excludes deselected data from
the first exported file. These are component and synthetic-workflow checks, not
general solver validation.
