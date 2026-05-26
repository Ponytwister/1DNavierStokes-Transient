#pragma once
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include <optimization.h>
#include <solvers.h>
#include <thread>
#include <stdafx.h>
#include <sqlite3.h>
#include <algorithm>
#include <map>
#include <assert.h>

// 7+ off
// 6: fatal
// 5: error
// 4: warn
// 3: info
// 2: debug
// 1: trace
// 0: all
struct report {
    int debug_level = 0;
    std::string report_text;
    report(int debug_level, std::string report_text) 
        : debug_level(debug_level)
        , report_text(std::move(report_text))
    {}
};

struct parameter_solution_struct {
    std::string PARAMETER;
    std::string VALUE;
    std::string UNITS;
    parameter_solution_struct(std::string PARAMETER, std::string VALUE, std::string UNITS) 
        : PARAMETER(std::move(PARAMETER))
        , VALUE(std::move(VALUE))
        , UNITS(std::move(UNITS)) 
    {}
};

struct solution_ID_struct {
    std::string SOLVE_SETTING_ID;
    std::string EXPERIMENT_NAME;
    std::string INLET_COND_ID;
    solution_ID_struct(int SOLVE_SETTING_ID, std::string EXPERIMENT_NAME, int INLET_COND_ID) 
        : SOLVE_SETTING_ID(std::to_string(SOLVE_SETTING_ID))
        , EXPERIMENT_NAME(std::move(EXPERIMENT_NAME))
        , INLET_COND_ID(std::to_string(INLET_COND_ID))
    {}
};

struct specie_struct {
    double diffusion_rate;
    //double difusion_dye; // m2/s From 4.9 × 10−6 cm2 s−1 The diffusion coefficient of fluorescein in water at 21.5°C, as calculated from the Wilke-Chang correlation
    //double difusion_beads; // m2/s From kB*T/(3*pi*visc*d) kB=1.380649×10−23 J⋅K−1
    double diameter;
    double particle_density;
    double molecular_weight;
    double QE;
    double r;
    int type; // 0 undefined, 1 molecule, 2 particle
    std::string input_units;
    std::string model_units;
    std::string name;
};

struct reaction_struct {
    double k[2];
    std::string name;
    std::map<specie_struct*, double> coef;
    std::map<specie_struct*, double> exp;
    std::vector<specie_struct*> specie_vect;
};

struct entrance_struct {
    std::map<specie_struct*, double> CONC;
    double ENTRANCE_FLOWRATE;
};

struct experiment_run_struct {
    double total_flowrate = 0.0d;
    double dye_conc_mgml;                  //= 0.00336d; // mg/ml FITC
    double dye_conc;                       //= dye_conc_mgml * 1000.0d / 332.326d * 6.022e+23;     // molecules FITC / m3
    double left_edge;
    double right_edge;
    double dt; // seconds
    double visc = 0.0010016d; // Dynamic viscosity of water at 20C in Pa.s
    double temperature = 20.0d + 273.15d;
    double solution_density = 1.0d;
    int number_of_reactions;
    int number_of_species;
    int low_ref_start, low_ref_end, high_ref_start, high_ref_end;
    ptrdiff_t FITC;
    ptrdiff_t PS_beads;
    ptrdiff_t Bound_Dye_1;
    ptrdiff_t Bound_Dye_2;
    ptrdiff_t FITC_Bead_1;
    ptrdiff_t FITC_Bead_2;
    std::string normalization_method;
    std::string name;
    std::vector<double> ENTRANCE_FLOWRATE;
    std::vector<specie_struct> species;
    std::vector<reaction_struct> reactions;
    std::vector<std::string> solve_for;
};

struct experiment_struct {
    double*species_in;
    double**species_out;
    int beginning_of_channel;
    int end_of_channel;
    int window_size;
    int window_start;
    int SOLUTION_ID = 0;
    int INLET_COND_ID = 0;
    double exp_DA;
    double model_DA;
    double exp_integral;
    double model_integral;
    double scale_factor;
    std::string second_name;
    std::vector<entrance_struct> entrances; 
    std::vector<double> model_profile;
    std::vector<double> experimental_profile;
    std::vector<double> error;
    std::vector<double> numeric_derivative;
    std::vector<double> experimental_derivative;
    experiment_run_struct* run;
    experiment_struct* zero_row_ptr;
    experiment_run_struct* zero_run;
};
void set_inlet_conc(experiment_struct* exp_ptr, double* solution);

typedef struct parameters_struct {
    // Device Dimenssions
    const double W = 5e-4, H = 4e-5, L = 0.025;  //meters: 500 um, 40 um, 2.5 cm
    int debug_level = 0;
    std::vector<report> state; // text output
    // Model Control Parameters
    int SOLVE_SETTING_ID = 0;
    bool SOLVE_SETTING_RECURSIVE_CALL = false;
    bool SOLUTION_ID_RECURSIVE_CALL = false;
    bool PROFILE_IDs_RECURSIVE_CALL = false;
    int Z, X, exp_left_padding, exp_right_padding;
    alglib::ae_int_t max_iterations;
    int iterations = 0;
    std::string output_file_name, scatter_correction_type;
    bool run_solver = true, disable_reactions = false, save_model_profiles = true, disable_reverse_reactions = false, use_alglib_init_values = true;
    double convergence_epsx;
    double time_step_convergence = 0.00001d;

    // Experimental Parameters
    int number_of_variables = 0;
    std::vector<std::string> solve_for;
    std::vector<std::string> global_solve_for;
    std::vector<std::string> retrieved;
    std::vector<experiment_run_struct> experiment_runs;
    std::vector<experiment_struct> experiments;
    int row_count = 0;
    int total_window_size = 0;

    // Alglib Inputs
    double* initial_values_alglib;
    double* low_bound;
    double* up_bound; 
    double* scale;
} parameters_t;

// Text FILE HANDLERS
void save_excel_output(std::string file_name);

// MODEL
void normalize_profile();
void model(const alglib::real_1d_array &control_parameters, alglib::real_1d_array &residuals, int row);
void alglib_solver(const alglib::real_1d_array &control_parameters, alglib::real_1d_array &residuals, void *ptr);
double scattering_correction(double NS, double species_2, double p, experiment_run_struct* run_ptr);
double& variable_location(const std::string& variable_name, experiment_run_struct* run_ptr);
void clear_output(int debug_level, std::string text_to_clear);
void add_report(int debug_level, std::string text_to_add);
void pop_report(int debug_level);
void pop_and_add(int debug_level, std::string text_to_add);
void finish_report(int debug_level, std::string text_to_add);
void add_finishing_report(int debug_level, std::string text_to_add);
void removeSpaces(std::string &str);
ptrdiff_t specie_index(experiment_run_struct* run_ptr, std::string specie_to_find, std::string second_specie_name = "na");
ptrdiff_t reaction_index(experiment_run_struct* run_ptr, std::string reaction_to_find, std::string second_reaction_name = "na");
double unit_conversion(experiment_run_struct* run_ptr, int specie, std::string current_unit, std::string desired_unit);
double lin_interpolate(double x, double x1, double y1, double x2, double y2);
std::string double_to_string(double arg);

// DB HANDLERS
void lines_from_profile_text(sqlite3* db);
void delete_values_from_db(sqlite3* db, std::string table, std::string where_conditions);
int exp_parameters_db_callback(void *data, int count, char **argv, char **columnNames);
void read_exp_parameters_from_db(sqlite3* db);
int read_model_parameters_db_callback(void *data, int count, char **argv, char **columnNames);
void read_model_parameters_from_db(sqlite3* db);
int raw_profiles_db_callback(void *data, int count, char **argv, char **columnNames);
void read_raw_profiles_from_db(sqlite3* db);
int inlet_cond_db_callback(void *data, int count, char **argv, char **columnNames);
void read_inlet_cond_from_db(sqlite3* db);
int get_SOLUTION_ID_from_db_callback(void *data, int count, char **argv, char **columnNames);
void get_SOLUTION_IDs_from_db(sqlite3* db);
int alglib_input_db_callback(void *data, int count, char **argv, char **columnNames);
void read_alglib_values_from_db(sqlite3* db);
void get_solve_settings_ID_from_db(sqlite3* db);
void write_model_profile_to_db(sqlite3* db);
void write_alglib_values_to_db(sqlite3* db);
int specie_db_callback(void *data, int count, char **argv, char **columnNames);
int reaction_db_callback(void *data, int count, char **argv, char **columnNames);
void read_specie_and_reaction_values_from_db(sqlite3* db);

extern parameters_t p;

struct model_profile_struct {
    std::string SOLUTION_ID;
    std::string X;
    std::string Free_Dye;
    std::string Bound_Dye;
    std::string Total_Dye;
    std::string Unbound_Beads;
    std::string Bound_Beads;
    std::string Total_Beads;
    std::string Experimental_Derivative;
    std::string Numeric_Derivative;
    std::string Experimental;
    std::string Numeric;
    model_profile_struct(int SOLUTION_ID, double X) 
        : SOLUTION_ID(std::to_string(SOLUTION_ID))
        , X(double_to_string(X)) 
    {}
};