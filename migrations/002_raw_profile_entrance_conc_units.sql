-- Per-profile units for ENTRANCE_CONC. NULL means use the FITC entry in the
-- experiment's SPECIE_MODEL_CONC_UNITS list.
BEGIN IMMEDIATE;
ALTER TABLE raw_profile ADD COLUMN ENTRANCE_CONC_UNITS TEXT;
COMMIT;
