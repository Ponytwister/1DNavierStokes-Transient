-- Allow raw profile inlet concentrations to store #variable_name aliases.
-- Numeric values remain numeric. Existing values are copied unchanged.
BEGIN IMMEDIATE;
CREATE TABLE raw_profile_entrance_concentrations_v5 (
    NAME TEXT NOT NULL,
    WT_PERCENT TEXT NOT NULL,
    ENTRANCE_NUMBER INTEGER NOT NULL
        CHECK(typeof(ENTRANCE_NUMBER)='integer' AND ENTRANCE_NUMBER >= 1),
    SPECIES_NAME TEXT NOT NULL,
    CONCENTRATION ANY NOT NULL
        CHECK(
            (typeof(CONCENTRATION) IN ('real','integer') AND CONCENTRATION >= 0 AND CONCENTRATION <= 1.7976931348623157e308)
            OR
            (typeof(CONCENTRATION)='text' AND instr(CONCENTRATION, '#') > 0
                AND length(CONCENTRATION) > instr(CONCENTRATION, '#')
                AND instr(CONCENTRATION, ' ') = 0
                AND instr(CONCENTRATION, char(9)) = 0
                AND instr(CONCENTRATION, char(10)) = 0
                AND instr(CONCENTRATION, char(13)) = 0)
        ),
    UNITS TEXT,
    PRIMARY KEY(NAME, WT_PERCENT, ENTRANCE_NUMBER, SPECIES_NAME)
);
INSERT INTO raw_profile_entrance_concentrations_v5(
    NAME, WT_PERCENT, ENTRANCE_NUMBER, SPECIES_NAME, CONCENTRATION, UNITS)
SELECT NAME, WT_PERCENT, ENTRANCE_NUMBER, SPECIES_NAME, CONCENTRATION, UNITS
FROM raw_profile_entrance_concentrations;
DROP TABLE raw_profile_entrance_concentrations;
ALTER TABLE raw_profile_entrance_concentrations_v5 RENAME TO raw_profile_entrance_concentrations;
COMMIT;
