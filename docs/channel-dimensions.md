# Per-experiment channel dimensions

The Experiments tab exposes **Channel Width (m)**, **Channel Height (m)**, and
**Channel Length (m)** beside NAME. Use **Add...** to create an experiment or
select a row and choose **Modify...** to edit its metadata and dimensions. All
three dimensions must be finite positive numbers in meters. These are physical
channel dimensions, distinct from the existing `WIDTH` profile/fit parameter.
New experiments start with dimensions 5e-4 x 4e-5 x 0.025 m.

The dialog's **Save** writes the experiment in one transaction; **Cancel**
discards edits. Names cannot be changed when modifying an experiment. A save
rejects concurrent changes to the original row; cancel and **Refresh** to retry.
Saving invalidates old results. The table itself is read-only.
A setup JSON still stores paths and model controls only; channel dimensions come
from its referenced database at run time.

## Migration and recovery

Close other writers and back up the database. With Python 3.7 or newer (standard
library only), run from the repository root:

```powershell
python tools/migrate_channel_dimensions.py navier.db
```

An optional `--backup PATH` chooses the backup destination, which must not exist.
The script creates a consistent SQLite backup, then applies the versioned
`migrations/001_channel_dimensions.sql` transaction. All existing rows receive
`5e-4`, `4e-5`, and `0.025`, respectively. New rows receive the same defaults.
Constraints reject NULL, nonnumeric, nonpositive, and infinite values. The script
checks database integrity afterward. Repeating it preserves already migrated
values; a partial schema is rejected. The SQL file itself is a one-time migration.
No existing experiment, profile, fit, or result rows are replaced.

On failure, the migration transaction rolls back. To recover the previous schema,
close all connections and restore the printed backup file to the original path.
Keep any newer database separately if it contains later work. Do not restore a
backup over an open database. The migration does not repurpose `PRAGMA user_version`.

Unmigrated databases remain readable and use the historical dimensions; the tab
explains that migration is needed for editing. A partial or invalid geometry is
rejected before model solution records are created. Old application versions
ignore these new columns and continue using hardcoded dimensions: use the updated
application whenever dimensions differ from the defaults.

## Model behavior and validation

`channel_dimensions` owns `const double W, H, L`, initialized in its constructor.
`parameters_struct` retains these immutable legacy defaults through inheritance.
Each `experiment_run_struct` also inherits an immutable geometry initialized from
its database row before any profile or parameter links are created. A new session
reads saved dimensions again; editing cannot change an active calculation.

Combined runs may use different geometry for each selected experiment. Each uses
its own residence time `W*H*L / total_flowrate`, axial time step `residence_time/Z`,
diffusion coefficient `D*dt*X*X/(W*W)`, profile scaling, integration window, and
saved profile coordinates. Existing boundary conditions, units, and parameter
links are unchanged. Reports repeat the coordinate header when channel width
changes, so mixed-width profiles have correctly labeled coordinates in micrometers.

Regression cases compare default and doubled dimensions: at flow `5e-10 m^3/s`
and four axial steps, time steps are 0.25 s and 2 s. A constant inlet with zero
reaction and no-flux walls remains uniform for both geometries, providing an
independent solution. Tests retain finite-output checks and the existing tolerance
`1e-12 + 1e-10*abs(expected)`. Compile-time assertions verify const dimensions.
These component and synthetic workflow cases are not validation of every transient
or experimental regime, and no performance claim is made.

Run the additional migration/backup tests with:

```powershell
python tests/channel_migration_test.py
```
