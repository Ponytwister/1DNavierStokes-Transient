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
