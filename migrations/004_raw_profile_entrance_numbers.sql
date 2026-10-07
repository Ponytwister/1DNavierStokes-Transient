-- Track a concentration override by both species and inlet.
-- Apply migration 003 first; existing per-species rows belong to inlet 1.
BEGIN IMMEDIATE;
CREATE TABLE raw_profile_entrance_concentrations_v4 (
    NAME TEXT NOT NULL,
    WT_PERCENT TEXT NOT NULL,
    ENTRANCE_NUMBER INTEGER NOT NULL
        CHECK(typeof(ENTRANCE_NUMBER)='integer' AND ENTRANCE_NUMBER >= 1),
    SPECIES_NAME TEXT NOT NULL,
    CONCENTRATION REAL NOT NULL
        CHECK(typeof(CONCENTRATION) IN ('real','integer') AND CONCENTRATION >= 0 AND CONCENTRATION <= 1.7976931348623157e308),
    UNITS TEXT,
    PRIMARY KEY(NAME, WT_PERCENT, ENTRANCE_NUMBER, SPECIES_NAME)
);
INSERT INTO raw_profile_entrance_concentrations_v4(
    NAME, WT_PERCENT, ENTRANCE_NUMBER, SPECIES_NAME, CONCENTRATION, UNITS)
SELECT NAME, WT_PERCENT, 1, SPECIES_NAME, CONCENTRATION, UNITS
FROM raw_profile_entrance_concentrations;
DROP TABLE raw_profile_entrance_concentrations;
ALTER TABLE raw_profile_entrance_concentrations_v4 RENAME TO raw_profile_entrance_concentrations;
COMMIT;
