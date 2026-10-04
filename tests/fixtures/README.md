# Independent regression fixtures

`inlet_profile.txt` first line contains: cell count, molecular weight (g/mol), inlet A
flow, inlet B flow (same arbitrary flow units), inlet A concentration, inlet B concentration
(mg/ml). Remaining values are the expected cell concentrations in the code's `umol`
unit convention (micromoles/liter for this conversion).

The four-cell case has equal inlet flows. A 100 g/mol solute at 0.002 mg/ml is 0.002 g/l,
or 20 micromoles/liter. The first two cells therefore contain 20 and the last two zero.
Both the current and previous solution buffers must be initialized to `[20, 20, 0, 0]`.
Expected values are hand-derived, not recorded from the implementation.

`model_controls.sql` creates only the two-column table needed by
`read_model_parameters_from_db`. Tests load it into SQLite `:memory:` and close the
connection afterwards. It is not a full application database or a copy of `navier.db`.

Floating-point comparisons require finite values and use absolute tolerance 1e-12 plus
relative tolerance 1e-10 times the magnitude of the expected value. These tolerances are
for the small helper/input cases here, not a proposed accuracy guarantee for full solves.
The interpolation case independently uses y = 3x + 2; the unit conversion case uses a
100 g/mol solute. Future full-solver cases should document their own physical assumptions,
boundary conditions, units, reference solution, and justified tolerances.

## Full workflow fixture

`workflow.sql` is a hand-authored, synthetic database schema and input set, not a
copy of `navier.db` or a production migration. Its column order follows the current
SQLite callback interface. Four profiles satisfy the inlet loader's existing
fourth-profile diagnostic. Their secondary labels identify rows only; no bead
concentration is supplied. There is no `0.0` reference row, so this case does not
test comparisons against another model worker's zero-profile output.

The channel is 500 micrometers wide, 40 micrometers high, and 0.025 meters long
(the model's existing dimensions), with eight cells and four time steps. A single
inlet supplies the same dye concentration to every cell. The configured volume
flow is 5e-10 m^3/s, giving a residence time of one second and a 0.25-second step.
At molecular weight 100 g/mol, 0.002 mg/ml = 0.002 g/l = 20 micromoles/liter
(the model's `umol` convention). Beads and bound dye start at zero. The configured
forward rate is zero, so both forward and reverse rates are zero. With no net
reaction and a spatially constant concentration, diffusion with the existing
no-flux edges preserves [20, 0, 0] for the three species at every cell and time.
The dye quantum efficiency is one, so the model profile is 20 and the normalized
experimental profile is one, giving zero residuals. These expectations are
analytical, not regenerated output. The unchanged finite-value tolerance is
`1e-12 + 1e-10 * abs(expected)`.

`keq1` is the single optimizer parameter, initially one in the reaction record.
It is intentionally unidentifiable when the forward rate is zero: this case
exercises real ALGLIB/model callbacks but does not demonstrate parameter recovery.
The separate `alglib_input` value 0.75 makes explicit fitted-input persistence
observable; the loader's established link/initial-value semantics select the
reaction value. Tests require a positive termination code and actual residual
evaluations, not a fixed platform-sensitive evaluation count.

Nine raw-profile samples include the channel's right endpoint, while the model
has eight cells. Saving therefore produces nine positions per profile (36 total);
the endpoint has experimental/model values but no additional concentration cell.
The text export's free-dye rows must each have eight values of 0.002 mg/ml.
Sentinel records with ID 99 belong to an unrelated experiment and must survive
saves. Repeated profile saves replace the selected results without duplication;
loading again reuses the four solution identities and solve-settings record.

Workflow tests also cancel after the first real model evaluation, confirm no
result profiles or fitted inputs were saved, and successfully run again. A zero
optimizer scale tests partial-load failure and recovery with corrected inputs.
The CLI fixture helper creates its own database and checks both save choices.
No production schema, numerical equation, tolerance, or dependency changes are
part of these workflow tests.

`ScatterControl.DatabaseAndInMemoryChoicesReachResiduals` adds uniform beads to
this fixture at the Gaussian center of each existing scatter calibration, with
particle diameters of 20 and 40 nm. Zero reaction rate and a uniform inlet keep
the dye profile at 20 and the bead concentration constant. The independently
simplified correction is `1 - amplitude + center * slope`; the expected squared
residual is `(correction - 1)^2`, or zero with `none`. The test covers database
controls and in-memory controls that deliberately contradict the database,
including an off/on/off sequence of fresh sessions. It uses the existing finite
value tolerance. The model profile itself remains uncorrected under both choices;
the correction currently affects only residuals. This checks control propagation
and the existing formula, not the empirical calibration or general transient
solver accuracy.
