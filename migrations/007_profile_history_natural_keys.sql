-- Replace inlet-ID based solve/result tables with append-only profile history.
-- Historical outputs and normalized inlet rows are intentionally discarded;
-- raw_profile itself is left untouched.
PRAGMA foreign_keys=OFF;
BEGIN IMMEDIATE;

DROP TABLE IF EXISTS model_profile;
DROP TABLE IF EXISTS parameter_solutions;
DROP TABLE IF EXISTS solutions;
DROP TABLE IF EXISTS solve_settings;
DROP TABLE IF EXISTS inlet_conditions;
DROP TABLE IF EXISTS profile_run_parameters;
DROP TABLE IF EXISTS profile_results;
DROP TABLE IF EXISTS profile_runs;

DELETE FROM raw_profile_entrance_concentrations;

CREATE TABLE profile_calculations (
    RUN_AT_UTC_US INTEGER NOT NULL PRIMARY KEY,
    ALL_EXP_FITTED TEXT,
    PARAMETERS_SOLVED_FOR TEXT,
    REACTIONS_ENABLED TEXT,
    SCATTER_METHOD TEXT,
    X_RESOLUTION INTEGER,
    Z_RESOLUTION INTEGER
);

CREATE TABLE profile_records (
    RUN_AT_UTC_US INTEGER NOT NULL,
    EXPERIMENT_NAME TEXT NOT NULL,
    WT_PERCENT TEXT NOT NULL,
    EXP_DA NUMERIC,
    MODEL_DA NUMERIC,
    EXP_INTEGRAL NUMERIC,
    MODEL_INTEGRAL NUMERIC,
    PRIMARY KEY (RUN_AT_UTC_US, EXPERIMENT_NAME, WT_PERCENT),
    FOREIGN KEY (RUN_AT_UTC_US) REFERENCES profile_calculations(RUN_AT_UTC_US) ON DELETE CASCADE
);
CREATE INDEX profile_records_profile_run
    ON profile_records(EXPERIMENT_NAME, WT_PERCENT, RUN_AT_UTC_US);

CREATE TABLE profile_record_parameters (
    RUN_AT_UTC_US INTEGER NOT NULL,
    EXPERIMENT_NAME TEXT NOT NULL DEFAULT '',
    WT_PERCENT TEXT NOT NULL DEFAULT '',
    PARAMETER TEXT NOT NULL,
    VALUE NUMERIC,
    UNITS TEXT,
    SOURCE TEXT NOT NULL DEFAULT '',
    PRIMARY KEY (RUN_AT_UTC_US, EXPERIMENT_NAME, WT_PERCENT, PARAMETER, SOURCE),
    FOREIGN KEY (RUN_AT_UTC_US) REFERENCES profile_calculations(RUN_AT_UTC_US) ON DELETE CASCADE
);

CREATE TABLE profile_record_samples (
    RUN_AT_UTC_US INTEGER NOT NULL,
    EXPERIMENT_NAME TEXT NOT NULL,
    WT_PERCENT TEXT NOT NULL,
    X NUMERIC NOT NULL,
    Free_Dye NUMERIC,
    Bound_Dye NUMERIC,
    Total_Dye NUMERIC,
    Unbound_Beads NUMERIC,
    Bound_Beads NUMERIC,
    Total_Beads NUMERIC,
    Experimental_Derivative NUMERIC,
    Numeric_Derivative NUMERIC,
    Experimental NUMERIC,
    "Numeric" NUMERIC,
    Error NUMERIC,
    Experimental_Difference NUMERIC,
    Numeric_Difference NUMERIC,
    PRIMARY KEY (RUN_AT_UTC_US, EXPERIMENT_NAME, WT_PERCENT, X),
    FOREIGN KEY (RUN_AT_UTC_US, EXPERIMENT_NAME, WT_PERCENT)
        REFERENCES profile_records(RUN_AT_UTC_US, EXPERIMENT_NAME, WT_PERCENT) ON DELETE CASCADE
);
CREATE INDEX profile_record_samples_lookup
    ON profile_record_samples(EXPERIMENT_NAME, WT_PERCENT, RUN_AT_UTC_US, X);

PRAGMA user_version=3;
COMMIT;
PRAGMA foreign_keys=ON;
