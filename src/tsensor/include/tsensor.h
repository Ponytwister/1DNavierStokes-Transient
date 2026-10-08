#pragma once
#include <channel_dimensions.h>
#include <cmath>
#include <string>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <optimization.h>
#include <solvers.h>
#include <thread>
#include <stdafx.h>
#include <sqlite3.h>
#include <algorithm>
#include <map>
#include <assert.h>
#include <deque>
#include <mutex>
#include <exception>
#include <feedback.h>
#include <stop_token>

enum level {
  all,
  trace,
  debug,
  info,
  warn,
  error,
  fatal,
  off
};
// 7+ off
// 6: fatal
// 5: error
// 4: warn
// 3: info
// 2: debug
// 1: trace
// 0: all

struct solvable {
    double default_value = 0.0;
    solvable* source = nullptr;
    bool param_init = false;
    std::string name;
    std::string source_name;

    solvable() = default;

    explicit solvable(double value)
        : default_value(value)
    {}

    solvable(double value, std::string name)
        : default_value(value)
        , name(std::move(name))
    {}

    solvable(std::string name)
        : name(std::move(name))
    {}

    solvable(std::string name, std::string source_name)
        : name(std::move(name))
        , source_name(std::move(source_name))
    {}

    explicit solvable(solvable* source_ptr, std::string name = "")
        : source(source_ptr)
        , name(std::move(name))
    {}

    double& value()
    {
        solvable* current = this;

        for (int i = 0; i < 100; ++i)
        {
            if (!current->source)
                return current->default_value;

            current = current->source;
        }

        throw std::runtime_error("Cycle detected in solvable chain");
    }

    bool is_linked() const
    {
        return source != nullptr;
    }

    void unlink()
    {
        if (source)
        {
            default_value = source->value();
            source = nullptr;
        }
    }

    void update()
    {
        if (source)
        {
            default_value = source->value();
        }
    }

    void init()
    {
        solvable* current = this;
        bool cycle_detected = true;
        for (int i = 0; i < 100; ++i)
        {
            current->param_init = true;
            if (!current->source) {
                cycle_detected = false;
                break;
            }
            current = current->source;
        }
        if (cycle_detected) {
            throw std::runtime_error("Cycle detected in solvable chain");
        }
    }
};

struct report {
    int debug_level = 0;
    std::string report_text;
    report(int debug_level, std::string report_text) 
        : debug_level(debug_level)
        , report_text(std::move(report_text))
    {}
};

struct specie_struct {
    double diffusion_rate;
    std::string diffusion_alias;
    //double difusion_dye; // m2/s From 4.9 × 10−6 cm2 s−1 The diffusion coefficient of fluorescein in water at 21.5°C, as calculated from the Wilke-Chang correlation
    //double difusion_beads; // m2/s From kB*T/(3*pi*visc*d) kB=1.380649×10−23 J⋅K−1
    double diameter;
    std::string diameter_alias;
    double particle_density;
    double molecular_weight;
    solvable QE;
    double r;
    int type; // 0 undefined, 1 molecule, 2 particle
    std::string input_units;
    std::string model_units;
    std::string name;
};

struct parameter_alias {
    std::string name;
    double sign = 1.0;
};

struct reaction_struct {
    solvable k[2];
    std::string k_alias[2];
    std::string name;
    std::map<specie_struct*, solvable> coef;
    std::map<specie_struct*, parameter_alias> coef_alias;
    std::map<specie_struct*, solvable> exp;
    std::map<specie_struct*, parameter_alias> exp_alias;
    std::vector<specie_struct*> specie_vect;
};

struct entrance_struct {
    std::map<specie_struct*, double> CONC;
    double ENTRANCE_FLOWRATE;
};

struct experiment_run_struct : channel_dimensions {
    explicit experiment_run_struct(channel_dimensions dimensions = channel_dimensions{})
        : channel_dimensions(dimensions) {}
    double total_flowrate = 0.0d;
    double dye_conc_mgml = 0;                  //= 0.00336d; // mg/ml FITC
    double dye_conc = 0;                       //= dye_conc_mgml * 1000.0d / 332.326d * 6.022e+23;     // molecules FITC / m3
    solvable left_edge;
    solvable width;
    double dt = 0; // seconds
    double visc = 0.0010016d; // Dynamic viscosity of water at 20C in Pa.s
    double temperature = 20.0d + 273.15d;
    double solution_density = 1.0d;
    int number_of_reactions = 0;
    int number_of_species = 0;
    int low_ref_start = 0, low_ref_end = 0, high_ref_start = 0, high_ref_end = 0;
    ptrdiff_t FITC = 0;
    ptrdiff_t PS_beads = 0;
    ptrdiff_t Bound_Dye_1 = 0;
    ptrdiff_t Bound_Dye_2 = 0;
    ptrdiff_t FITC_Bead_1 = 0;
    ptrdiff_t FITC_Bead_2 = 0;
    std::string normalization_method;
    std::string name;
    std::vector<double> ENTRANCE_FLOWRATE;
    std::vector<specie_struct> species;
    std::vector<reaction_struct> reactions;
    std::map<std::string, solvable> alias_variables;
    std::vector<std::string> solve_for;
};

struct experiment_struct {
    std::vector<std::vector<double>> species_out;
    int beginning_of_channel;
    int end_of_channel;
    int window_size;
    int window_start;
    int SOLUTION_ID = 0;
    int INLET_COND_ID = 0;
    bool has_legacy_inlet_cond_id = false;
    bool has_entrance_conc_override = false;
    double entrance_conc_override = 0.0;
    std::string entrance_conc_alias;
    std::string entrance_conc_units;
    struct entrance_concentration_override {
        int entrance_number;
        std::string species_name;
        double concentration;
        std::string units;
    };
    std::vector<entrance_concentration_override> entrance_concentrations;
    double exp_DA;
    double model_DA;
    double exp_integral;
    double model_integral;
    double analytic_exp_integral;
    double analytic_model_integral;
    double scale_factor;
    bool omit;
    solvable left_edge;
    solvable width;
    std::string second_name;
    std::vector<entrance_struct> entrances;
    std::vector<double> channel_position;
    std::vector<double> model_profile;
    std::vector<double> raw_experimental_profile;
    std::vector<double> experimental_profile;
    std::vector<double> error;
    std::vector<double> numeric_derivative;
    std::vector<double> experimental_derivative;
    std::vector<double> numeric_difference;
    std::vector<double> experimental_difference;
    std::vector<double> analytical_zero;
    experiment_run_struct* run = nullptr; // Non-owning, within this session.
    experiment_struct* zero_row_ptr = nullptr;
};
struct parameters_struct;
void set_inlet_conc(parameters_struct& p, experiment_struct* exp_ptr, double* solution);

typedef struct parameters_struct : channel_dimensions {
    explicit parameters_struct(channel_dimensions defaults = channel_dimensions{})
        : channel_dimensions(defaults) {}
    // Links point into this state's containers. Copying or moving would leave
    // pointers referring to the old owner; create a fresh state for each run.
    parameters_struct(const parameters_struct&) = delete;
    parameters_struct& operator=(const parameters_struct&) = delete;
    parameters_struct(parameters_struct&&) = delete;
    parameters_struct& operator=(parameters_struct&&) = delete;
    // Inherited const W/H/L are legacy defaults. Each experiment owns its geometry.
    int debug_level = 0;
    std::vector<report> state; // text output
    std::recursive_mutex report_mutex; // Shared by the existing model workers.
    tsensor_workflow::progress_callback progress;
    std::exception_ptr progress_failure;
    std::stop_token cancellation;
    tsensor_workflow::operation active_operation = tsensor_workflow::operation::none;
    // Model Control Parameters
    int SOLVE_SETTING_ID = 0;
    bool SOLVE_SETTING_RECURSIVE_CALL = false;
    bool SOLUTION_ID_RECURSIVE_CALL = false;
    bool PROFILE_IDs_RECURSIVE_CALL = false;
    int Z = 0, X = 0, exp_left_padding = 0, exp_right_padding = 0;
    alglib::ae_int_t max_iterations = 0;
    int iterations = 0;
    std::string output_file_name, scatter_correction_type = "none";
    bool run_solver = true, disable_reactions = false, save_model_profiles = true, disable_reverse_reactions = false, use_alglib_init_values = true;
    double convergence_epsx = 1e-9;
    double time_step_convergence = 0.0001d;

    // Experimental Parameters
    std::vector<std::string> solve_for;
    std::vector<std::string> global_solve_for;
    std::vector<std::string> retrieved;
    std::deque<solvable> solvables;
    std::vector<experiment_run_struct> experiment_runs;
    std::vector<experiment_struct> experiments;
    int row_count = 0;
    int total_window_size = 0;

    // Alglib Inputs
    std::vector<double> initial_values_alglib;
    std::vector<double> low_bound;
    std::vector<double> up_bound;
    std::vector<double> scale;

    std::map<std::string, double> initial_values_alglib_map;
    std::map<std::string, double> low_bound_map;
    std::map<std::string, double> up_bound_map; 
    std::map<std::string, double> scale_map;
} parameters_t;

// Callback deliveries are serialized and synchronous, possibly on a model worker.
// No callback means no console output. Do not reenter/mutate a running session.
void publish_event(parameters_t& p, tsensor_workflow::progress_event event);

// Text FILE HANDLERS
void save_excel_output(parameters_t& p, const std::filesystem::path& file_name);
// Extended, full-precision in-memory source for customizable desktop reports.
std::string generate_excel_report(parameters_t& p);

// MODEL
void normalize_profile(parameters_t& p);
void model(parameters_t& p, const alglib::real_1d_array &control_parameters, alglib::real_1d_array &residuals, int row);
void alglib_solver(const alglib::real_1d_array &control_parameters, alglib::real_1d_array &residuals, void *ptr);
double scattering_correction(parameters_t& p, double NS, double species_2, double coef, experiment_run_struct* run_ptr);
solvable& variable_location(const std::string& variable_name, experiment_run_struct* run_ptr);
void add_report(parameters_t& p, int debug_level, std::string text_to_add);
void pop_report(parameters_t& p, int debug_level);
void pop_and_add(parameters_t& p, int debug_level, std::string text_to_add);
void add_finishing_report(parameters_t& p, int debug_level, std::string text_to_add);
void removeSpaces(std::string &str);
ptrdiff_t specie_index(experiment_run_struct* run_ptr, std::string specie_to_find, std::string second_specie_name = "na");
ptrdiff_t reaction_index(experiment_run_struct* run_ptr, std::string reaction_to_find, std::string second_reaction_name = "na");
double unit_conversion(experiment_run_struct* run_ptr, int specie, std::string current_unit, std::string desired_unit);
double lin_interpolate(double x, double x1, double y1, double x2, double y2);
std::string double_to_string(double arg);

// DB HANDLERS
void lines_from_profile_text(parameters_t& p, sqlite3* db);
void delete_values_from_db(parameters_t& p, sqlite3* db, std::string table, std::string where_conditions);
int exp_parameters_db_callback(void *data, int count, char **argv, char **columnNames);
void read_exp_parameters_from_db(parameters_t& p, sqlite3* db);
int read_model_parameters_db_callback(void *data, int count, char **argv, char **columnNames);
void read_model_parameters_from_db(parameters_t& p, sqlite3* db);
int raw_profiles_db_callback(void *data, int count, char **argv, char **columnNames);
void read_raw_profiles_from_db(parameters_t& p, sqlite3* db);
int inlet_cond_db_callback(void *data, int count, char **argv, char **columnNames);
void read_inlet_cond_from_db(parameters_t& p, sqlite3* db);
int get_SOLUTION_ID_from_db_callback(void *data, int count, char **argv, char **columnNames);
void get_SOLUTION_IDs_from_db(parameters_t& p, sqlite3* db);
int alglib_input_db_callback(void *data, int count, char **argv, char **columnNames);
void read_alglib_values_from_db(parameters_t& p, sqlite3* db);
void get_solve_settings_ID_from_db(parameters_t& p, sqlite3* db);
int get_solvable_initial_values_from_db_callback(void *data, int count, char **argv, char **columnNames);
void get_solvable_initial_values_from_db(parameters_t& p, sqlite3* db);
void write_model_profile_to_db(parameters_t& p, sqlite3* db);
void write_alglib_values_to_db(parameters_t& p, sqlite3* db);
int specie_db_callback(void *data, int count, char **argv, char **columnNames);
int reaction_db_callback(void *data, int count, char **argv, char **columnNames);
void read_specie_and_reaction_values_from_db(parameters_t& p, sqlite3* db);

// Cooperative checkpoint; only the stop source may be used concurrently.
void check_cancellation(const parameters_t& p);

// Initialize per-experiment const dimensions before creating any model links.
void initialize_channel_dimensions(parameters_t& p, sqlite3* db);
