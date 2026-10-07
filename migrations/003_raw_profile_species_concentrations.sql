-- One concentration override per species and raw profile.
-- Apply migration 002 first. Existing scalar ENTRANCE_CONC values become FITC rows.
BEGIN IMMEDIATE;
CREATE TABLE raw_profile_entrance_concentrations (
    NAME TEXT NOT NULL,
    WT_PERCENT TEXT NOT NULL,
    SPECIES_NAME TEXT NOT NULL,
    CONCENTRATION REAL NOT NULL
        CHECK(typeof(CONCENTRATION) IN ('real','integer') AND CONCENTRATION >= 0 AND CONCENTRATION <= 1.7976931348623157e308),
    UNITS TEXT,
    PRIMARY KEY(NAME, WT_PERCENT, SPECIES_NAME)
);
INSERT INTO raw_profile_entrance_concentrations(NAME, WT_PERCENT, SPECIES_NAME, CONCENTRATION, UNITS)
SELECT NAME, WT_PERCENT, 'FITC', ENTRANCE_CONC, ENTRANCE_CONC_UNITS
FROM raw_profile
-- Older databases may store a multiline entrance concentration matrix in this
-- column. Migrate numeric scalar values only; keep other legacy values intact.
WHERE typeof(ENTRANCE_CONC) IN ('integer', 'real');
COMMIT;
