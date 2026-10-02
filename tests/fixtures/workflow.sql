-- Synthetic uniform equilibrium; no data from navier.db.
CREATE TABLE model_controls (criterion TEXT PRIMARY KEY, value TEXT);
INSERT INTO model_controls VALUES
 ('debug_level','7'), ('experiment_name','uniform'),
 ('width resolution (X)','8'), ('length/time resolution (Z)','4'),
 ('universal_solve_for','keq1'), ('max_iterations','3'),
 ('save_model_profiles','true'), ('scatter_correction_type','none');
CREATE TABLE experiments (NAME TEXT PRIMARY KEY, LOW_REF_LEFT INTEGER, LOW_REF_RIGHT INTEGER,
 HIGH_REF_LEFT INTEGER, HIGH_REF_RIGHT INTEGER, DEFAULT_NORMALIZATION TEXT,
 SPECIES TEXT, REACTIONS TEXT, ENTRANCE_FLOWRATE TEXT, EDGES TEXT,
 SPECIE_INLET_CONC_UNITS TEXT, SPECIE_MODEL_CONC_UNITS TEXT, WIDTH REAL);
INSERT INTO experiments VALUES ('uniform',0,0,1,8,'none',
 'FITC PS_40nm 40nm_Bound_Dye_1','FITC_40nm_1','0.0000000005','0',
 'mg/ml wt% mg/ml','umol umol umol',8);
CREATE TABLE species (SPECIES_NAME TEXT PRIMARY KEY, SPECIES_TYPE TEXT, DIFFUSION_RATE REAL,
 QE REAL, PARTICLE_DIAMETER REAL, PARTICLE_DENSITY REAL, MOLECULAR_WEIGHT REAL);
INSERT INTO species VALUES ('FITC','molecule',4.9e-10,1,0,1,100),
 ('PS_40nm','particle',0,0,40,1,100), ('40nm_Bound_Dye_1','molecule',4.9e-10,1,0,1,100);
CREATE TABLE reactions (REACTION_NAME TEXT PRIMARY KEY, SPECIES TEXT, COEFFICIENTS TEXT, Ks TEXT, EXPONENTS TEXT);
INSERT INTO reactions VALUES ('FITC_40nm_1','FITC PS_40nm 40nm_Bound_Dye_1','-1 -1 1','0 1','1 1 1');
CREATE TABLE raw_profile (NAME TEXT, WT_PERCENT TEXT, CHANNEL_LEFT_EDGE INTEGER,
 CHANNEL_RIGHT_EDGE INTEGER, INTENSITY_ARRAY TEXT, INLET_COND_ID INTEGER,
 LEFT_EDGE REAL, WIDTH REAL, OMIT TEXT);
WITH ids(id) AS (VALUES(1),(2),(3),(4))
INSERT INTO raw_profile SELECT 'uniform', printf('%d.0',id),0,9,
 '1'||char(9)||'1'||char(9)||'1'||char(9)||'1'||char(9)||'1'||char(9)||'1'||char(9)||'1'||char(9)||'1'||char(9)||'1',
 id,0,8,'false' FROM ids;
CREATE TABLE inlet_conditions (INLET_COND_ID INTEGER, SPECIE_CONC REAL, ENTRANCE_NUMBER INTEGER, SPECIES_NAME TEXT);
INSERT INTO inlet_conditions VALUES (1,0.002,1,'FITC'),(2,0.002,1,'FITC'),(3,0.002,1,'FITC'),(4,0.002,1,'FITC');
CREATE TABLE alglib_input (VARIABLE TEXT PRIMARY KEY, "INITIAL VALUE" REAL, "LOWER BOUND" REAL, "UPPER BOUND" REAL, SCALE REAL);
INSERT INTO alglib_input VALUES ('keq1',0.75,0.5,2,1);
CREATE TABLE solve_settings (SOLVE_SETTING_ID INTEGER PRIMARY KEY, ALL_EXP_FITTED TEXT,
 PARAMETERS_SOLVED_FOR TEXT, REACTIONS_ENABLED TEXT, SCATTER_METHOD TEXT, X_RESOLUTION INTEGER, Z_RESOLUTION INTEGER);
CREATE TABLE solutions (SOLUTION_ID INTEGER PRIMARY KEY, SOLVE_SETTING_ID INTEGER, EXPERIMENT_NAME TEXT,
 INLET_COND_ID INTEGER, EXP_DA REAL, MODEL_DA REAL, EXP_INTEGRAL REAL, MODEL_INTEGRAL REAL, SECOND_NAME TEXT);
CREATE TABLE parameter_solutions (SOLVE_SETTING_ID INTEGER, PARAMETER TEXT, VALUE REAL, UNITS TEXT, SOURCE TEXT, SOLUTION_ID INTEGER);
CREATE TABLE model_profile (SOLUTION_ID INTEGER, X REAL, Free_Dye REAL, Bound_Dye REAL, Total_Dye REAL,
 Unbound_Beads REAL, Bound_Beads REAL, Total_Beads REAL, Experimental_Derivative REAL,
 Numeric_Derivative REAL, Experimental REAL, Numeric REAL, Error REAL,
 Experimental_Difference REAL, Numeric_Difference REAL);
-- Unrelated records must survive every save.
INSERT INTO solve_settings VALUES (99,'unrelated','keq1','all reactions enabled','none',8,4);
INSERT INTO solutions VALUES (99,99,'unrelated',99,3,4,5,6,'untouched');
INSERT INTO parameter_solutions VALUES (99,'keq1',42,'','global',99);
INSERT INTO model_profile (SOLUTION_ID,X,Numeric) VALUES (99,0,42);
