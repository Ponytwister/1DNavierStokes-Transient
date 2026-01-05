#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include <eigen-3.4.0/Eigen/Dense>
#include <optimization.h>
#include <stdafx.h>
#include <solvers.h>
#include <sqlite3.h>
#include <thread>
// D: Dye
// NS: Nanoparticle Site
// species_2: Nanoparticle Bound Dye

using namespace std;
// Text FILE HANDLERS
void prime_excel_output(string file_name);
void save_excel_output();

// DB HANDLERS
void lines_from_profile_text(sqlite3* db);
void delete_values_from_db(sqlite3* db, string table, string where_conditions);
int exp_parameters_db_callback(void *data, int count, char **argv, char **columnNames);
void read_exp_parameters_from_db(sqlite3* db);
int read_model_parameters_db_callback(void *data, int count, char **argv, char **columnNames);
void read_model_parameters_from_db(sqlite3* db);
int raw_profiles_db_callback(void *data, int count, char **argv, char **columnNames);
void read_raw_profiles_from_db(sqlite3* db);
int alglib_input_db_callback(void *data, int count, char **argv, char **columnNames);
void read_alglib_values_from_db(sqlite3* db);
void write_normalized_values_to_db(sqlite3* db);
void write_model_profile_to_db(sqlite3* db);
void write_alglib_values_to_db(sqlite3* db);
int specie_db_callback(void *data, int count, char **argv, char **columnNames);
int reaction_db_callback(void *data, int count, char **argv, char **columnNames);
void read_specie_and_reaction_values_from_db(sqlite3* db);

// MODEL
void normalize_profile();
void model(const alglib::real_1d_array &control_parameters, alglib::real_1d_array &residuals, int row);
void alglib_solver(const alglib::real_1d_array &control_parameters, alglib::real_1d_array &residuals, void *ptr);
double scattering_correction(double NS, double species_2, double p);
double& variable_location(const std::string& variable_name);

struct parameters_struct {
    // Device Dimenssions
    const double W = 5e-4, H = 4e-5, L = 0.025;  //meters: 500 um, 40 um, 2.5 cm 

    // Model Control Parameters
    int Z, X, exp_left_padding, exp_right_padding;
    alglib::ae_int_t max_iterations;
    int iterations = 0;
    string experiment_name, file_name, output_file_name, scatter_correction_type;;
    bool run_solver = true, disable_reactions = false, save_normalized_profiles = true, save_model_profiles = true;
    double convergence_epsx;

    // Experimental Parameters
    string normalization_method;
    double*ENTRANCE_FLOWRATE;
    double visc = 0.0010016d, temperature = 20.0d + 273.15d; // Dynamic viscosity of water at 20C in Pa.s
    double diameter_old;                    // meters
    double* wt_percent;                     //= 0.1d; // wt%
    double* dye_conc_mgml;                  //= 0.00336d; // mg/ml FITC
    double* dye_conc;                       //= dye_conc_mgml * 1000.0d / 332.326d * 6.022e+23;     // molecules FITC / m3
    int low_ref_start, low_ref_end, high_ref_start, high_ref_end;
    int number_of_variables, number_of_species, number_of_reactions, number_of_entrances;
    std::vector<std::string> solve_for, reaction, SPECIE;
    bool solve_for_left_pad = false, solve_for_right_pad = false;
    bool solve_for_p1 = false, solve_for_kon1 = false, solve_for_koff1 = false, solve_for_keq1 = false, solve_for_QE1 = false;
    bool solve_for_p2 = false, solve_for_kon2 = false, solve_for_koff2 = false, solve_for_keq2 = false, solve_for_QE2 = false;
    double left_pad, right_pad, p1, kon1, koff1, keq1, QE1, p2, kon2, koff2, keq2, QE2;

    // Species 
    std::vector<std::string> specie_type;
    double* diffusion_rate;
    double* diameter;
    double* QE;
    double difusion_dye; // m2/s From 4.9 × 10−6 cm2 s−1 The diffusion coefficient of fluorescein in water at 21.5°C, as calculated from the Wilke-Chang correlation
    double difusion_beads; // m2/s From kB*T/(3*pi*visc*d) kB=1.380649×10−23 J⋅K−1
    double* species_unit_conversion;        // mass / umolar conversion, FITC mgml / umolar conversion, wt% / umolar conversion
    double dt; // seconds
    double* r;
    Eigen::Matrix<double, -1, -1>** inverted_diffusion_matrix;

    // Reactions
    std::vector<std::string> react_specie_text;
    std::vector<std::string> react_coef_text;
    double** coef;
    std::vector<std::string> react_ks_text;
    double* kon;
    double* koff;
    double* keq;
    std::vector<std::string> react_exp_text;
    double** exp;

    // Eperimental Profile
    double***CONC;
    double**species_in;
    double**species_out;
    int* beginning_of_channel;
    int* end_of_channel;
    double* experimental_profile;
    int row_count = 0;
    int window_size = 180;
    double num_profile_width;
    double* numeric_model_profile;

    // Alglib Inputs
    double* initial_values_alglib;
    double* low_bound;
    double* up_bound; 
    double* scale;
} parameters;

void
prime_excel_output(string file_name)
{
    ofstream fout;
    fout.open(file_name, std::ofstream::out | std::ofstream::trunc);
    fout << "res_time  "    << "bind_ratio(p1)  "   << "forward_reaction_rate_1  "  << "reverse_reaction_rate_1  "  << "dye_conc.   "   << "bead_conc.  "   << "QuenchEnhance_Factor_1  "   << "species  "; 
    for (int x = 0; x < parameters.X; x++) {
        fout << x << "  ";
    }
    fout << endl;
    fout << "sec  "         << "molc/bead  "        << "forward_reaction_rate_1  "  << "reverse_reaction_rate_1  "  << "molc/m3   "     << "bead/m3  "      << "ND/D_Ratio  "               << "Channel_Width_(um)->  "; 
}

void
save_excel_output()
{
    double dye_bead_ratio = 2000.0d;
    double total_dye = 0.0d;
    //int z = parameters.Z - 1;
    ofstream fout;
    fout.open(parameters.output_file_name, std::ofstream::out | std::ofstream::app);

    double scale_factor = parameters.W * 1.0e6 / ((double)parameters.window_size - parameters.left_pad - parameters.right_pad);
    for (int x = 0; x < parameters.X; x++) {
        fout << (x  - parameters.left_pad) * scale_factor << "  ";
    }
    fout << endl;

    double num_derivative[parameters.window_size];
    double exp_derivative[parameters.window_size];
    num_derivative[0] = 0;
    exp_derivative[0] = 0;
    num_derivative[parameters.window_size-1] = 0;
    exp_derivative[parameters.window_size-1] = 0;
    
    for (int row = 0; row < parameters.row_count; row++) {
        for (int i = 1; i < parameters.window_size-1; i++) {
            num_derivative[i] = (parameters.numeric_model_profile[row * parameters.window_size + i + 1] / parameters.dye_conc[row] - parameters.numeric_model_profile[row * parameters.window_size + i - 1] / parameters.dye_conc[row]) / (2 * parameters.W * 1.0e+6 / parameters.window_size);
            exp_derivative[i] = (parameters.experimental_profile[row * parameters.window_size + i + 1] - parameters.experimental_profile[row * parameters.window_size + i - 1]) / (2 * parameters.W * 1.0e+6 / parameters.window_size);
        }
        for (int j = 0; j < 10; j++) {
            fout << parameters.Z * parameters.dt << "  " << parameters.p1 << "  " << parameters.kon1 << "  " << parameters.koff1 << "  " << parameters.dye_conc[row] << "   " <<  parameters.wt_percent[row] << "  " <<  parameters.QE1 << "  ";
            switch(j) {
                case 0:
                    fout << "Free_Dye  ";
                    for (int x = 0; x < parameters.X; x++) {
                        fout << parameters.species_out[0][x + row * parameters.X] / parameters.dye_conc[row] << "    ";
                    }
                    break;
                case 1:
                    fout << "Bound_Dye  ";
                    for (int x = 0; x < parameters.X; x++) {
                        fout << (parameters.species_out[2][x + row * parameters.X] + parameters.species_out[3][x + row * parameters.X]) / parameters.dye_conc[row] << "    ";
                    }
                    break;
                case 2:
                    fout << "Total_Dye  ";
                    for (int x = 0; x < parameters.X; x++) {
                        fout << (parameters.species_out[0][x + row * parameters.X] + parameters.species_out[2][x + row * parameters.X] + parameters.species_out[3][x + row * parameters.X]) / parameters.dye_conc[row] << "    ";
                    }
                    break;
                case 3:
                    fout << "Unbound_Beads  ";
                    for (int x = 0; x < parameters.X; x++) {
                        fout << parameters.species_out[1][x + row * parameters.X] * parameters.species_unit_conversion[1] << "    ";
                    }
                    break;
                case 4:
                    fout << "Bound_Dye  ";
                    for (int x = 0; x < parameters.X; x++) {
                        fout << (parameters.species_out[2][x + row * parameters.X] / parameters.p1 + parameters.species_out[3][x + row * parameters.X] / parameters.p1 / parameters.p2) * parameters.species_unit_conversion[1] << "    ";
                    }
                    break;   
                case 5:
                    fout << "Total_Beads  ";
                    for (int x = 0; x < parameters.X; x++) {
                        fout << (parameters.species_out[1][x + row * parameters.X] + parameters.species_out[2][x + row * parameters.X] / parameters.p1 + parameters.species_out[3][x + row * parameters.X] / parameters.p1 / parameters.p2) * parameters.species_unit_conversion[1] << "    ";
                    }
                    break;
                case 6:
                    fout << "Experimental_Derivative  ";
                    for (int x = 0; x < parameters.window_size; x++) {
                        fout << exp_derivative[x] << "    ";
                    }
                    break;
                case 7:
                    fout << "Numeric_Derivative  ";
                    for (int x = 0; x < parameters.window_size; x++) {
                        fout << num_derivative[x] << "    ";
                    }
                    break;    
                case 8:
                    fout << "Experimental_Profile  ";
                    for (int x = 0; x < parameters.window_size; x++) {
                        fout << parameters.experimental_profile[row * parameters.window_size + x] << "    ";
                    }
                    break;
                case 9:
                    fout << "Total_Dye_rescale  ";
                    for (int x = 0; x < parameters.window_size; x++) {
                        fout << parameters.numeric_model_profile[row * parameters.window_size + x] / parameters.dye_conc[row] << "    ";
                    }
                    break;
            }
            fout << endl;
        }
        fout << endl;
    }
    fout.close(); 
}

int 
raw_profile_row_count_db_callback(void *data, int count, char **argv, char **columnNames)
{
    //1st parameter of this function is received from 4th parameter of sqlite3_exec
    //count -> is the number of columns
    //columnNames ->  array of pointers to strings where each entry represents the name of corresponding result column as obtained
    //argv -> array of pointers to strings obtained as if from [sqlite3_column_text()]
    string criterion;
    for(int i = 0; i < count; i++) {
        criterion = columnNames[i];
        if (criterion == "NAME") {
            if (argv[i] != NULL) {
                parameters.row_count++;
            } else {}
        }
    };
    return 0;
}

void
lines_from_profile_text(sqlite3* db)
{
    char* errMsg = 0;
    string sqltext = "SELECT NAME FROM 'raw_profile' WHERE NAME='" + parameters.experiment_name + "';";
    const char* sql = sqltext.c_str();
    int rc = sqlite3_exec(db, sql, raw_profile_row_count_db_callback, 0, &errMsg);
    if(rc != SQLITE_OK){
        printf("Error in executing SQL: %s \n", errMsg);
        sqlite3_free(errMsg);
    }
}

void 
delete_values_from_db(sqlite3* db, string table, string where_conditions) //reading data using callback functions
{
    int row = 0;
    char* errMsg = 0;
    string sqltext;
    if(where_conditions.empty()) {
        return;
    }
    sqltext = sqltext + "DELETE FROM " + table + " WHERE " + where_conditions + ";";
    int rc = sqlite3_exec(db, sqltext.c_str(), 0, 0, &errMsg);
    if(rc != SQLITE_OK){
        printf("Error in executing deletion SQL: %s \n", errMsg);
        cout << sqltext << endl;
        sqlite3_free(errMsg);
    }
}

int 
read_model_parameters_db_callback(void *data, int count, char **argv, char **columnNames)
{
    //1st parameter of this function is received from 4th parameter of sqlite3_exec
    //count -> is the number of columns
    //columnNames ->  array of pointers to strings where each entry represents the name of corresponding result column as obtained
    //argv -> array of pointers to strings obtained as if from [sqlite3_column_text()]
    string criterion;
    string value = argv[1];
    criterion = argv[0];
    if (criterion == "width resolution (X)") {
        if (!value.empty()) {parameters.X = stoi(value);} else {parameters.X = 500;}
    } else if (criterion == "length/time resolution (Z)") {
        if (!value.empty()) {parameters.Z = stoi(value);} else {parameters.Z = 100;}
    } else if (criterion == "experiment_name") {
        if (!value.empty()) {parameters.experiment_name = value;} else {parameters.experiment_name = "24_3";}
    } else if (criterion == "exp_left_padding") {
        if (!value.empty()) {parameters.exp_left_padding = stoi(value);} else {parameters.exp_left_padding = 20;}
    } else if (criterion == "exp_right_padding") {
        if (!value.empty()) {parameters.exp_right_padding = stoi(value);} else {parameters.exp_right_padding = 20;}
    } else if (criterion == "disable_reactions") {
        if (!value.empty()) {if(value == "true") {parameters.disable_reactions = true;} else {parameters.disable_reactions = false;}
        } else {parameters.disable_reactions = false;}
    } else if (criterion == "run_solver") {
        if (!value.empty()) {if(value == "true") {parameters.run_solver = true;} else {parameters.run_solver = false;}
        } else {parameters.run_solver = true;}
    } else if (criterion == "scatter_correction_type") {
        if (!value.empty()) {parameters.scatter_correction_type = value;} else {parameters.scatter_correction_type = "none";}
    } else if (criterion == "save_normalized_profiles") {
        if (!value.empty()) {if(value == "true") {parameters.save_normalized_profiles = true;} else {parameters.save_normalized_profiles = false;}
        } else {parameters.save_normalized_profiles = true;}
    } else if (criterion == "save_model_profiles") {
        if (!value.empty()) {if(value == "true") {parameters.save_model_profiles = true;} else {parameters.save_model_profiles = false;}
        } else {parameters.save_model_profiles = true;}
    } else if (criterion == "max_iterations") {
        if (!value.empty()) {parameters.max_iterations = stoi(value);} else {parameters.max_iterations = 0;}
    } else if (criterion == "convergence_epsx") {
        if (!value.empty()) {parameters.convergence_epsx = stod(value);} else {parameters.convergence_epsx = 1e-9;}
    }
    return 0;
}

void 
read_model_parameters_from_db(sqlite3* db) //reading data using callback functions
{
    char* errMsg = 0;
    string sqltext = "SELECT * FROM 'model_controls';";
    const char* sql = sqltext.c_str();
    int rc = sqlite3_exec(db, sql, read_model_parameters_db_callback, 0, &errMsg);
    if(rc != SQLITE_OK){
        printf("Error in executing SQL: %s \n", errMsg);
        sqlite3_free(errMsg);
    }
}

int 
exp_parameters_db_callback(void *data, int count, char **argv, char **columnNames)
{
    //1st parameter of this function is received from 4th parameter of sqlite3_exec
    //count -> is the number of columns
    //columnNames ->  array of pointers to strings where each entry represents the name of corresponding result column as obtained
    //argv -> array of pointers to strings obtained as if from [sqlite3_column_text()]
    string criterion;
    
    for(int i = 0; i < count; i++) {
        criterion = columnNames[i];
        if (criterion == "DYE_CONC_MGML") {
            for (int row = 0; row < parameters.row_count; row++) {
                if (argv[i] != NULL) {
                    parameters.dye_conc_mgml[row] = stod(argv[i]);
                    parameters.dye_conc[row] = parameters.dye_conc_mgml[row] / 332.326d * 1.0e+6;  // um mol     * 1000.0d / 332.326d * 6.022e+23;  // molecules / m3  
                } else {
                    parameters.dye_conc_mgml[row] = 0.0d;
                    parameters.dye_conc[row] = 0.0d;
                }
            }
        } else if (criterion == "BEAD_DIAMETER_NM") {
            if (argv[i] != NULL) {
                parameters.diameter_old = stod(argv[i]) * 1.0e-9d;
            } else {
                parameters.diameter_old = 41.0e-9d;
            }
        } else if (criterion == "LOW_REF_LEFT") {
            if (argv[i] != NULL) {
                parameters.low_ref_start = stoi(argv[i]);
            } else {
                parameters.low_ref_start = 0;
            }
        } else if (criterion == "LOW_REF_RIGHT") {
            if (argv[i] != NULL) {
                parameters.low_ref_end = stoi(argv[i]);
            } else {
                parameters.low_ref_end = 10;
            }
        } else if (criterion == "HIGH_REF_LEFT") {
            if (argv[i] != NULL) {
                parameters.high_ref_start = stoi(argv[i]);
            } else {
                parameters.high_ref_start = 95;
            }
        } else if (criterion == "HIGH_REF_RIGHT") {
            if (argv[i] != NULL) {
                parameters.high_ref_end = stoi(argv[i]);
            } else {
                parameters.high_ref_end = 105;
            }
        } else if (criterion == "PARAMETERS_TO_SOLVE_FOR") {
            string PARAMETERS_TO_SOLVE_FOR = argv[i];
            string s; // variable to store token obtained from the original string
            stringstream ss(PARAMETERS_TO_SOLVE_FOR); // constructing stream from the string
            while (getline(ss, s, ' ')) { 
                parameters.solve_for.push_back(s); // store token string in the vector
            }
            parameters.number_of_variables = parameters.solve_for.size();
            for (int item = 0; item < parameters.number_of_variables; item++) {
                if (argv[i] != NULL) {
                    if (parameters.solve_for[item] == "left_pad") {
                        parameters.solve_for_left_pad = true;
                    } else if (parameters.solve_for[item] == "right_pad") {
                        parameters.solve_for_right_pad = true;
                    } else if (parameters.solve_for[item] == "p1") {
                        parameters.solve_for_p1 = true;
                    } else if (parameters.solve_for[item] == "kon1") {
                        parameters.solve_for_kon1 = true;
                    } else if (parameters.solve_for[item] == "koff1") {
                        parameters.solve_for_koff1 = true;
                    } else if (parameters.solve_for[item] == "keq1") {
                        parameters.solve_for_keq1 = true;
                    } else if (parameters.solve_for[item] == "QE1") {
                        parameters.solve_for_QE1 = true;
                    } else if (parameters.solve_for[item] == "p2") {
                        parameters.solve_for_p2 = true;
                    } else if (parameters.solve_for[item] == "kon2") {
                        parameters.solve_for_kon2 = true;
                    } else if (parameters.solve_for[item] == "koff2") {
                        parameters.solve_for_koff2 = true;
                    } else if (parameters.solve_for[item] == "keq2") {
                        parameters.solve_for_keq2 = true;
                    } else if (parameters.solve_for[item] == "QE2") {
                        parameters.solve_for_QE2 = true;
                    }
                } else {

                }
            }
        } else if (criterion == "DEFAULT_NORMALIZATION") {
            if (argv[i] != NULL) {
                parameters.normalization_method = argv[i];
            } else {
                parameters.normalization_method = "none";
            }
        } else if (criterion == "SPECIES") {
            string SPECIES = argv[i];
            string s; // variable to store token obtained from the original string
            stringstream ss(SPECIES); // constructing stream from the string
            while (getline(ss, s, ' ')) { 
                parameters.SPECIE.push_back(s); // store token string in the vector
            }
            parameters.number_of_species = parameters.SPECIE.size();
        } else if (criterion == "REACTIONS") {
            string reactions = argv[i];
            string s; // variable to store token obtained from the original string
            stringstream ss(reactions); // constructing stream from the string
            while (getline(ss, s, ' ')) { 
                parameters.reaction.push_back(s); // store token string in the vector
            }
            parameters.number_of_reactions = parameters.reaction.size();
        } else if (criterion == "ENTRANCE_FLOWRATE") {
            string entrances = argv[i];
            string s; // variable to store token obtained from the original string
            stringstream ss(entrances); // constructing stream from the string
            std::vector<string> entrance_vector;
            while (getline(ss, s, ' ')) { 
                entrance_vector.push_back(s); // store token string in the vector
            }
            parameters.number_of_entrances = entrance_vector.size();
            parameters.ENTRANCE_FLOWRATE = new double[parameters.number_of_entrances];
            for (int entrance = 0; entrance < parameters.number_of_entrances; entrance++) {
                parameters.ENTRANCE_FLOWRATE[entrance] = stod(entrance_vector[entrance]);
            }
        }
        
    }
    return 0;
}

void 
read_exp_parameters_from_db(sqlite3* db) //reading data using callback functions
{
    char* errMsg = 0;
    string sqltext = "SELECT * FROM 'experiments' WHERE NAME='" + parameters.experiment_name + "';";
    const char* sql = sqltext.c_str();
    int rc = sqlite3_exec(db, sql, exp_parameters_db_callback, 0, &errMsg);
    if(rc != SQLITE_OK){
        printf("Error in executing SQL: %s \n", errMsg);
        sqlite3_free(errMsg);
    }
}

int 
raw_profiles_db_callback(void *data, int count, char **argv, char **columnNames)
{
    //1st parameter of this function is received from 4th parameter of sqlite3_exec
    //count -> is the number of columns
    //columnNames ->  array of pointers to strings where each entry represents the name of corresponding result column as obtained
    //argv -> array of pointers to strings obtained as if from [sqlite3_column_text()]
    string criterion;
    double wt_percentages[parameters.row_count];
    for (int test_row = 0; test_row < parameters.row_count; test_row++) {
        wt_percentages[test_row] = parameters.wt_percent[test_row];
    }
    double current_wt_percent;
    int col_index;
    int row = 0;
    double last_profile_point = 0;
    for(int i = 0; i < count; i++) {
        criterion = columnNames[i];
        if (criterion == "WT_PERCENT") {
            if (argv[i] != NULL) {
                current_wt_percent = stod(argv[i]);
                for (row = 0 ;row < parameters.row_count; row++) {
                    if ((wt_percentages[row] == 0) && (current_wt_percent == 0.0 || row != 0)) { 
                        break;
                    }
                }
                parameters.wt_percent[row] = current_wt_percent;
            } else {}
        } else if (criterion == "CHANNEL_LEFT_EDGE") {
            if (argv[i] != NULL) {
                parameters.beginning_of_channel[row] = stoi(argv[i]);
            } else {parameters.beginning_of_channel[row] = 57;}
        } else if (criterion == "CHANNEL_RIGHT_EDGE") {
            if (argv[i] != NULL) {
                parameters.end_of_channel[row] = stoi(argv[i]);
            } else {parameters.end_of_channel[row] = 168;}
        } else if (criterion == "INTENSITY_ARRAY") {
            string INTENSITY_ARRAY = argv[i];
            string s; // variable to store token obtained from the original string
            stringstream ss(INTENSITY_ARRAY); // constructing stream from the string
            vector<string> intensity; // declaring vector to store the string after split
            // using while loop until the getline condition is satisfied
            // ' ' represent split the string whenever a space is found in the original string 
            while (getline(ss, s, '	')) { 
                intensity.push_back(s); // store token string in the vector
            }
            for (int x = 0; x < intensity.size(); x++) {
                col_index = x - parameters.beginning_of_channel[row] + parameters.exp_left_padding;
                if (col_index < 180) {
                    if (argv[i] != NULL) {
                        if (col_index < 0) { 
                            // Ignore these data points
                        } else if (col_index < parameters.exp_left_padding) { // Points added in left padding
                            parameters.experimental_profile[row * parameters.window_size + col_index] = 0.0;
                        } else if (col_index < parameters.end_of_channel[row] + parameters.beginning_of_channel[row] - parameters.exp_left_padding) { // Actual data 
                            parameters.experimental_profile[row * parameters.window_size + col_index] = stod(intensity[x]);
                            last_profile_point = stod(intensity[x]);
                        } else if (col_index < parameters.end_of_channel[row] + parameters.beginning_of_channel[row] - parameters.exp_left_padding) { // Points added in right padding
                            parameters.experimental_profile[row * parameters.window_size + col_index] = last_profile_point;
                        }
                    } else {parameters.experimental_profile[row * parameters.window_size + col_index] = 1;}
                }
            }
        } else if (criterion == "ENTRANCE_CONC") {
            if (argv[i] != NULL) {
                string entrances = argv[i];
                string s; // variable to store token obtained from the original string
                stringstream ss(entrances); // constructing stream from the string
                string s2;
                int entrance = 0;
                while (getline(ss, s)) {
                    stringstream ss2(s);
                    int specie = 0;
                    while (getline(ss2, s2, ' ')) {
                        parameters.CONC[entrance][specie][row] = stod(s2);
                        specie = specie + 1;
                    }
                    entrance = entrance + 1;
                }
            } else {
                throw std::runtime_error("No Entrance Concentrations");}
        }
    };
    return 0;
}

void 
read_raw_profiles_from_db(sqlite3* db) //reading data using callback functions
{
    char* errMsg = 0;
    string sqltext = "SELECT * FROM 'raw_profile' WHERE NAME='" + parameters.experiment_name + "';";
    const char* sql = sqltext.c_str();
    int rc = sqlite3_exec(db, sql, raw_profiles_db_callback, 0, &errMsg);
    if(rc != SQLITE_OK){
        printf("Error in executing SQL: %s \n", errMsg);
        sqlite3_free(errMsg);
    }
}

int 
alglib_input_db_callback(void *data, int count, char **argv, char **columnNames)
{
    //1st parameter of this function is received from 4th parameter of sqlite3_exec
    //count -> is the number of columns
    //columnNames ->  array of pointers to strings where each entry represents the name of corresponding result column as obtained
    //argv -> array of pointers to strings obtained as if from [sqlite3_column_text()]
    string criterion;
    string current_variable;
    bool varible_to_solve_for;
    int index = 0;
    for(int i = 0; i < count; i++) {
        criterion = columnNames[i];
        if (criterion == "VARIABLE") {
            if (argv[i] != NULL) {
                current_variable = argv[i];
                //varible_to_solve_for = (parameters.solve_for_left_pad == true && current_variable == "left_pad") || (parameters.solve_for_right_pad == true && current_variable == "right_pad" ) || (parameters.solve_for_p1 == true && current_variable == "p1") || (parameters.solve_for_kon1 == true && current_variable == "kon1") || (parameters.solve_for_koff1 == true && current_variable == "koff1") || (parameters.solve_for_keq1 == true && current_variable == "keq1") || (parameters.solve_for_QE1 == true && current_variable == "QE1") || (parameters.solve_for_p2 == true && current_variable == "p2") || (parameters.solve_for_kon2 == true && current_variable == "kon2") || (parameters.solve_for_koff2 == true && current_variable == "koff2") || (parameters.solve_for_keq2 == true && current_variable == "keq2") || (parameters.solve_for_QE2 == true && current_variable == "QE2");
                varible_to_solve_for = (find(parameters.solve_for.begin(), parameters.solve_for.end(), current_variable) != parameters.solve_for.end());
                for (int item = 0; item < parameters.number_of_variables; item++) {
                    if (parameters.solve_for[item] == current_variable) {
                        index = item;
                    }
                }
            } else {
                throw std::runtime_error("No variable name");
            }
        } else if (criterion == "INITIAL VALUE") {
            if (argv[i] != NULL) {
                variable_location(current_variable) = stod(argv[i]);
                if (varible_to_solve_for) {
                    parameters.initial_values_alglib[index] = stod(argv[i]); // initial values for alglib solver
                };
            } else {
                throw std::runtime_error("No Initial Value for" + current_variable);
            }
        } else if (criterion == "LOWER BOUND") {
            if (argv[i] != NULL) {
                if (varible_to_solve_for) {
                    parameters.low_bound[index] = stod(argv[i]);
                }
            } else {
                throw std::runtime_error("No Lower Bound for" + current_variable);
            }
        } else if (criterion == "UPPER BOUND") {
            if (argv[i] != NULL) {
                if (varible_to_solve_for) {
                    parameters.up_bound[index] = stod(argv[i]);
                }
            } else {
                throw std::runtime_error("No Upper Bound for" + current_variable);
            }
        } else if (criterion == "SCALE") {
            if (argv[i] != NULL) {
                if (varible_to_solve_for) {
                    parameters.scale[index] = stod(argv[i]);
                }
            } else {
                throw std::runtime_error("No Scale for" + current_variable);
            }
        } 
    };
    return 0;
}

void 
read_alglib_values_from_db(sqlite3* db) //reading data using callback functions
{
    char* errMsg = 0;
    string sqltext = "SELECT * FROM 'alglib_input';";
    int rc = sqlite3_exec(db, sqltext.c_str(), alglib_input_db_callback, 0, &errMsg);
    if(rc != SQLITE_OK){
        printf("Error in executing SQL: %s \n", errMsg);
        sqlite3_free(errMsg);
    }
}

int 
specie_db_callback(void *data, int count, char **argv, char **columnNames)
{
    //1st parameter of this function is received from 4th parameter of sqlite3_exec
    //count -> is the number of columns
    //columnNames ->  array of pointers to strings where each entry represents the name of corresponding result column as obtained
    //argv -> array of pointers to strings obtained as if from [sqlite3_column_text()]
    string criterion;
    ptrdiff_t specie_index;
    for(int i = 0; i < count; i++) {
        criterion = columnNames[i];
        if (criterion == "SPECIES_NAME") {
            if (argv[i] != NULL) {
                specie_index = distance(parameters.SPECIE.begin(), find(parameters.SPECIE.begin(), parameters.SPECIE.end(), argv[i]));
            }
        } else if (criterion == "SPECIES_TYPE") {
            if (argv[i] != NULL) {
                parameters.specie_type[specie_index] = argv[i];
            } else {
                parameters.specie_type[specie_index] = "particle";
            }
        } else if (criterion == "DIFFUSION_RATE") {
            if (argv[i] != NULL) {
                parameters.diffusion_rate[specie_index] = stod(argv[i]);
            } else {
                parameters.diffusion_rate[specie_index] = 0.0d;
            }
        } else if (criterion == "QE") {
            if (argv[i] != NULL) {
                parameters.QE[specie_index] = stod(argv[i]);
            } else {
                parameters.QE[specie_index] = 0.0d;
            }
        } else if (criterion == "PARTICLE_DIAMETER") {
            if (argv[i] != NULL) {
                parameters.diameter[specie_index] = stod(argv[i]) * 1.0e-9d;
            } else {
                parameters.diameter[specie_index] = 0.0d;
            }
        }
    };
    return 0;
}

int 
reaction_db_callback(void *data, int count, char **argv, char **columnNames)
{
    //1st parameter of this function is received from 4th parameter of sqlite3_exec
    //count -> is the number of columns
    //columnNames ->  array of pointers to strings where each entry represents the name of corresponding result column as obtained
    //argv -> array of pointers to strings obtained as if from [sqlite3_column_text()]
    string criterion;
    ptrdiff_t reaction_index;
    std::vector<std::ptrdiff_t> specie_vect;
    std::vector<double> coef_vect;
    std::vector<double> ks_vect;
    std::vector<double> exp_vect;
    for(int i = 0; i < count; i++) {
        criterion = columnNames[i];
        if (criterion == "REACTION_NAME") {
            if (argv[i] != NULL) {
                reaction_index = distance(parameters.reaction.begin(), find(parameters.reaction.begin(), parameters.reaction.end(), argv[i]));
            }
        } else if (criterion == "SPECIES") {
            string species = argv[i];
            parameters.react_specie_text.at(reaction_index) = species;
            string s; // variable to store token obtained from the original string
            stringstream ss(species); // constructing stream from the string
            while (getline(ss, s, ' ')) { 
                specie_vect.push_back(distance(parameters.SPECIE.begin(), find(parameters.SPECIE.begin(), parameters.SPECIE.end(), s)));
            }
        } else if (criterion == "COEFFICIENTS") {
            string coefficients = argv[i];
            parameters.react_coef_text.at(reaction_index) = coefficients;
            string s; // variable to store token obtained from the original string
            stringstream ss(coefficients); // constructing stream from the string
            while (getline(ss, s, ' ')) { 
                coef_vect.push_back(stod(s)); // store token string in the vector
            }
        } else if (criterion == "Ks") {
            string rate_constants = argv[i];
            parameters.react_ks_text.at(reaction_index) =  rate_constants;
            string s; // variable to store token obtained from the original string
            stringstream ss(rate_constants); // constructing stream from the string
            while (getline(ss, s, ' ')) { 
                ks_vect.push_back(stod(s)); // store token string in the vector
            }
        } else if (criterion == "EXPONENTS") {
            string exponents = argv[i];
            parameters.react_exp_text.at(reaction_index) = exponents;
            string s; // variable to store token obtained from the original string
            stringstream ss(exponents); // constructing stream from the string
            while (getline(ss, s, ' ')) { 
                exp_vect.push_back(stod(s)); // store token string in the vector
            }
        }

    };
    for (int specie = 0; specie < parameters.number_of_species; specie++) {
        parameters.coef[reaction_index][specie] = 0;
    }
    for (int reactant = 0; reactant < specie_vect.size(); reactant++) {
        parameters.coef[reaction_index][specie_vect[reactant]] = coef_vect.at(reactant);
        parameters.exp[reaction_index][specie_vect[reactant]] = exp_vect.at(reactant);
    }
    parameters.kon[reaction_index] = ks_vect.at(0);
    parameters.keq[reaction_index] = ks_vect.at(1);
    return 0;
}

void 
read_specie_and_reaction_values_from_db(sqlite3* db) //reading data using callback functions
{
    char* errMsg = 0;
    string sqltext = "SELECT * FROM 'species' WHERE ";
    for (int item = 0; item < parameters.number_of_species; item++) {
        sqltext = sqltext + "SPECIES_NAME = '" + parameters.SPECIE[item] + "'";
        if (item != parameters.number_of_species - 1) {
            sqltext = sqltext + " OR ";
        } else {
            sqltext = sqltext + ";";
        }
    }
    int rc = sqlite3_exec(db, sqltext.c_str(), specie_db_callback, 0, &errMsg);
    if(rc != SQLITE_OK){
        printf("Error in executing SQL: %s \n", errMsg);
        sqlite3_free(errMsg);
    }
    sqltext = "SELECT * FROM 'reactions' WHERE ";
    for (int item = 0; item < parameters.number_of_reactions; item++) {
        sqltext = sqltext + "REACTION_NAME = '" + parameters.reaction[item] + "'";
        if (item != parameters.number_of_reactions - 1) {
            sqltext = sqltext + " OR ";
        } else {
            sqltext = sqltext + ";";
        }
    }
    rc = sqlite3_exec(db, sqltext.c_str(), reaction_db_callback, 0, &errMsg);
    if(rc != SQLITE_OK){
        printf("Error in executing SQL: %s \n", errMsg);
        sqlite3_free(errMsg);
    }
    
}

void 
write_normalized_values_to_db(sqlite3* db) //reading data using callback functions
{
    int row = 0;
    char* errMsg = 0;
    string sqltext;
    delete_values_from_db(db, "normalized_profile", "NAME = '" + parameters.experiment_name + "' AND NORM_METHOD = '" + parameters.normalization_method + "'");

    for (int row = 0; row < parameters.row_count; row++) {
        sqltext = sqltext + "INSERT INTO normalized_profile (NAME, WT_PERCENT, NORM_METHOD, INTENSITY_ARRAY) ";
        sqltext = sqltext + "VALUES ('" + parameters.experiment_name + "','" + std::to_string(parameters.wt_percent[row]) + "','" + parameters.normalization_method + "','";

        // INTENSITY_ARRAY
        for (int i = 0; i < 180; i++) {
            sqltext = sqltext + std::to_string(parameters.experimental_profile[row * parameters.window_size + i]) + " ";
        }
        sqltext = sqltext + "'); ";
    }
    int rc = sqlite3_exec(db, sqltext.c_str(), 0, 0, &errMsg);
    if(rc != SQLITE_OK){
        printf("Error in executing SQL: %s \n", errMsg);
        sqlite3_free(errMsg);
    }
}

void 
write_model_profile_to_db(sqlite3* db) //reading data using callback functions
{
    char* errMsg = 0;
    string sqltext;
    string species;
    string INTENSITY_ARRAY;
    int width = parameters.X;
    string LATERAL_POSITION_ARRAY;
    delete_values_from_db(db, "model_profile", "NAME = '" + parameters.experiment_name + "' AND NORM_METHOD = '" + parameters.normalization_method + "' AND SCATTER_METHOD = '" + parameters.scatter_correction_type + "'");
    double scale_factor = parameters.W * 1.0e6 / ((double)parameters.window_size - parameters.left_pad - parameters.right_pad);
    for (int x = 0; x < width; x++) {
        LATERAL_POSITION_ARRAY = LATERAL_POSITION_ARRAY + std::to_string((x  - parameters.left_pad) * scale_factor) + " ";
    }
    for (int row = 0; row < parameters.row_count; row++) {
        for (int specie = 0; specie < 12; specie++) {
            INTENSITY_ARRAY = "";
            switch(specie) {
                case 0:
                    species = "Free Dye";
                    width = parameters.X;
                    for (int x = 0; x < width; x++) {
                        INTENSITY_ARRAY = INTENSITY_ARRAY + std::to_string(parameters.species_out[0][x + row * parameters.X] / parameters.dye_conc[row]) + " ";
                    }
                    break;
                case 1:
                    species = "Bound Dye";
                    width = parameters.X;
                    for (int x = 0; x < width; x++) {
                        INTENSITY_ARRAY = INTENSITY_ARRAY + std::to_string(parameters.species_out[2][x + row * parameters.X] / parameters.dye_conc[row]) + " ";
                    }
                    break;
                case 2:
                    species = "Double Bound Dye";
                    width = parameters.X;
                    for (int x = 0; x < width; x++) {
                        INTENSITY_ARRAY = INTENSITY_ARRAY + std::to_string( (1 + 1 / parameters.p2) * parameters.species_out[3][x + row * parameters.X] / parameters.dye_conc[row]) + " ";
                    }
                    break;
                case 3:
                    species = "Total Dye";
                    width = parameters.X;
                    for (int x = 0; x < width; x++) {
                        INTENSITY_ARRAY = INTENSITY_ARRAY + std::to_string((parameters.species_out[0][x + row * parameters.X] + parameters.species_out[2][x + row * parameters.X] + (1 + 1 / parameters.p2) * parameters.species_out[3][x + row * parameters.X]) / parameters.dye_conc[row]) + " ";
                    }
                    break;
                case 4:
                    species = "Unbound Beads (wt%)";
                    width = parameters.X;
                    for (int x = 0; x < parameters.X; x++) {
                        INTENSITY_ARRAY = INTENSITY_ARRAY + std::to_string(parameters.species_out[1][x + row * parameters.X] * parameters.species_unit_conversion[1]) + " ";
                    }
                    break;
                case 5:
                    species = "Bound Beads (wt%)";
                    width = parameters.X;
                    for (int x = 0; x < width; x++) {
                        INTENSITY_ARRAY = INTENSITY_ARRAY + std::to_string(parameters.species_out[2][x + row * parameters.X] / parameters.p1 * parameters.species_unit_conversion[1]) + " ";
                    }
                    break;
                case 6:
                    species = "Double Bound Beads (wt%)";
                    width = parameters.X;
                    for (int x = 0; x < width; x++) {
                        INTENSITY_ARRAY = INTENSITY_ARRAY + std::to_string(parameters.species_out[3][x + row * parameters.X] / parameters.p1 / parameters.p2 * parameters.species_unit_conversion[1]) + " ";
                    }
                    break; 
                case 7:
                    species = "Total Beads (wt%)";
                    width = parameters.X;
                    for (int x = 0; x < width; x++) {
                        INTENSITY_ARRAY = INTENSITY_ARRAY + std::to_string((parameters.species_out[1][x + row * parameters.X] + parameters.species_out[2][x + row * parameters.X] / parameters.p1 + parameters.species_out[3][x + row * parameters.X] / parameters.p1 / parameters.p2) * parameters.species_unit_conversion[1]) + " ";
                    }
                    break;
                case 8:
                    species = "Experimental Derivative";
                    width = parameters.window_size;
                    for (int x = 0; x < (width - 1); x++) {
                        if ((x == 0) || (x + 1) == (width - 1)) {
                            INTENSITY_ARRAY = INTENSITY_ARRAY + "0 ";
                        } else {
                            INTENSITY_ARRAY = INTENSITY_ARRAY + std::to_string((parameters.experimental_profile[row * parameters.window_size + x + 1] - parameters.experimental_profile[row * parameters.window_size + x - 1]) / (2 * parameters.W * 1.0e+6 / parameters.window_size)) + " ";
                        }
                    }
                    break;
                case 9:
                    species = "Numeric Derivative";
                    width = parameters.window_size;
                    for (int x = 0; x < (width - 1); x++) {
                        if ((x == 0) || (x + 1) == (width - 1)) {
                            INTENSITY_ARRAY = INTENSITY_ARRAY + "0 ";
                        } else {
                            INTENSITY_ARRAY = INTENSITY_ARRAY + std::to_string((parameters.numeric_model_profile[row * parameters.window_size + x + 1] / parameters.dye_conc[row] - parameters.numeric_model_profile[row * parameters.window_size + x - 1] / parameters.dye_conc[row]) / (2 * parameters.W * 1.0e+6 / parameters.window_size)) + " ";
                        }
                    }
                    break;    
                case 10:
                    species = "Experimental Profile";
                    width = parameters.window_size;
                    for (int x = 0; x < width; x++) {
                        INTENSITY_ARRAY = INTENSITY_ARRAY + std::to_string(parameters.experimental_profile[row * parameters.window_size + x]) + " ";
                    }
                    break;
                case 11:
                    species = "Numeric Model Profile";
                    width = parameters.window_size;
                    for (int x = 0; x < width; x++) {
                        INTENSITY_ARRAY = INTENSITY_ARRAY + std::to_string(parameters.numeric_model_profile[row * parameters.window_size + x] / parameters.dye_conc[row]) + " ";
                    }
                    break;
            }
            sqltext = sqltext + "INSERT INTO model_profile (NAME, WT_PERCENT, NORM_METHOD, SCATTER_METHOD, SPECIES, left_pad, right_pad, p1, kon1, koff1, QE1, p2, kon2, koff2, QE2, keq1, keq2, LATERAL_POSITION_ARRAY, INTENSITY_ARRAY)";
            sqltext = sqltext + " VALUES ('" + parameters.experiment_name + "','" + std::to_string(parameters.wt_percent[row]) + "','" + parameters.normalization_method + "','" + parameters.scatter_correction_type + "','" + species + "','" + std::to_string(parameters.left_pad) + "','" + std::to_string(parameters.right_pad) + "','" + std::to_string(parameters.p1) + "','" + std::to_string(parameters.kon1) + "','" + std::to_string(parameters.koff1) + "','" + std::to_string(parameters.QE1) + "','" + std::to_string(parameters.p2) + "','" + std::to_string(parameters.kon2) + "','" + std::to_string(parameters.koff2) + "','" + std::to_string(parameters.QE2) + "','" + std::to_string(parameters.keq1) + "','" + std::to_string(parameters.keq2) + "','" + LATERAL_POSITION_ARRAY + "','" + INTENSITY_ARRAY + "'); ";
        }
        
    }

    int rc = sqlite3_exec(db, sqltext.c_str(), 0, 0, &errMsg);
    if(rc != SQLITE_OK){
        printf("Error in executing SQL: %s \n", errMsg);
        sqlite3_free(errMsg);
    }
}

void 
write_alglib_values_to_db(sqlite3* db) //reading data using callback functions
{
    char* errMsg = 0;
    string sqltext;
    for (int index = 0; index < parameters.number_of_variables; index++) {
        sqltext = sqltext + "UPDATE alglib_input SET 'INITIAL VALUE' = '" + std::to_string(parameters.initial_values_alglib[index]) + "' WHERE VARIABLE = '" + parameters.solve_for[index] + "'; ";
    }
    for (int specie = 0; specie < parameters.number_of_species; specie++) { 
        sqltext = sqltext + "UPDATE species SET 'QE' = '" + std::to_string(parameters.QE[specie]) + "' WHERE SPECIES_NAME = '" + parameters.SPECIE.at(specie) + "'; ";
    }
    for (int reaction = 0; reaction < parameters.number_of_reactions; reaction++) {
        ptrdiff_t reaction_index = distance(parameters.reaction.begin(), find(parameters.reaction.begin(), parameters.reaction.end(), parameters.reaction.at(reaction)));
        sqltext = sqltext + "UPDATE reactions SET 'Ks' = '" + std::to_string(parameters.kon[reaction_index]) + " " + std::to_string(parameters.keq[reaction_index]) + "' WHERE REACTION_NAME = '" + parameters.reaction.at(reaction_index) + "'; ";

        std::vector<std::ptrdiff_t> specie_vect;
        string s_react; // variable to store token obtained from the original string
        stringstream ss_react(parameters.react_specie_text.at(reaction_index)); // constructing stream from the string
        while (getline(ss_react, s_react, ' ')) { 
            specie_vect.push_back(distance(parameters.SPECIE.begin(), find(parameters.SPECIE.begin(), parameters.SPECIE.end(), s_react)));
        }
        sqltext = sqltext + "UPDATE reactions SET 'COEFFICIENTS' = '";
        for (std::ptrdiff_t specie = 0; specie < specie_vect.size(); specie++) {
            std::cout << parameters.coef[reaction_index][specie_vect.at(specie)];
            if (specie != specie_vect.size() - 1) {
                std::cout << " ";
            }
        }
        sqltext = sqltext + "' REACTION_NAME = '" + parameters.reaction.at(reaction_index) + "'; ";

        sqltext = sqltext + "UPDATE reactions SET 'EXPONENTS' = '";
        for (std::ptrdiff_t specie = 0; specie < specie_vect.size(); specie++) {
            std::cout << parameters.exp[reaction_index][specie_vect.at(specie)];
            if (specie != specie_vect.size() - 1) {
                std::cout << " ";
            }
        }
        sqltext = sqltext + "' REACTION_NAME = '" + parameters.reaction.at(reaction_index) + "'; ";
    }
    int rc = sqlite3_exec(db, sqltext.c_str(), 0, 0, &errMsg);
    if(rc != SQLITE_OK){
        printf("Error in executing SQL: %s \n", errMsg);
        sqlite3_free(errMsg);
    }
}

void 
normalize_profile() 
{
    double low_ref, high_ref, peak_ref;
    for (int j = 0; j < parameters.row_count; j++) {
        low_ref = 0.0d;
        high_ref = 0.0d;
        peak_ref = 0.0d;
        for (int i = parameters.low_ref_start + parameters.exp_left_padding; i < parameters.high_ref_end + parameters.exp_left_padding; i++) {
            if (i <= parameters.low_ref_end + parameters.exp_left_padding) {
                low_ref += parameters.experimental_profile[j * parameters.window_size + i];
            } else if (i >= parameters.high_ref_start + parameters.exp_left_padding) {
                high_ref += parameters.experimental_profile[j * parameters.window_size + i];
            }
        }
        for (int i = parameters.exp_left_padding; i < parameters.high_ref_end + parameters.exp_left_padding; i++) {
            peak_ref = max(parameters.experimental_profile[j * parameters.window_size + i], peak_ref);
        }
        if (parameters.low_ref_end < parameters.low_ref_start) {
            low_ref = 0;
        } else {
            low_ref = low_ref / (double)(parameters.low_ref_end - parameters.low_ref_start + 1);
        }
        
        high_ref = high_ref / (double)(parameters.high_ref_end - parameters.high_ref_start);
        for (int i = 0; i <= parameters.window_size; i++) {
            if (parameters.normalization_method == "ridge linear scaling") {
                if(i < parameters.low_ref_start + parameters.exp_left_padding) {
                    parameters.experimental_profile[j * parameters.window_size + i] = 0;
                } else if (i >= parameters.high_ref_end + parameters.exp_left_padding) {  
                    parameters.experimental_profile[j * parameters.window_size + i] = 1;
                } else {
                    parameters.experimental_profile[j * parameters.window_size + i] = (parameters.experimental_profile[j * parameters.window_size + i] - low_ref) / (high_ref - low_ref);
                }
            } else if (parameters.normalization_method == "peak linear scaling") {
                if (i < parameters.low_ref_start + parameters.exp_left_padding) {
                    parameters.experimental_profile[j * parameters.window_size + i] = 0;
                } else if(i >= parameters.window_size - parameters.exp_right_padding) {
                    parameters.experimental_profile[j * parameters.window_size + i] = (high_ref - low_ref) / (peak_ref - low_ref);
                } else {
                    parameters.experimental_profile[j * parameters.window_size + i] = (parameters.experimental_profile[j * parameters.window_size + i] - low_ref) / (peak_ref - low_ref);
                }
            } else if (parameters.normalization_method == "none") {

            } 
        }
    }
}

double
scattering_correction(double species_1, double species_2, double p)
{
    if (parameters.scatter_correction_type == "NS_ND") {
        double bead_wt = (species_1 + species_2 / p) * parameters.species_unit_conversion[1];
        if (parameters.diameter_old < 30.0e-9d) {
            return -0.238826108843563 * std::exp(-pow(0.0258645848310996 - bead_wt, 2.0d) / (2.0d * pow(0.00418992180425235, 2.0d))) + bead_wt * 1.40044738887211 + 1;
        };
        return -0.448191328804794 * std::exp(-pow(0.0711222856018783 - bead_wt, 2.0d) / (2.0d * pow(0.012623365952763, 2.0d))) + bead_wt * 2.14776259822044 + 1;
    } else {
        return 1.00d;
    }
}

double&
variable_location(const std::string& variable_name) {
    if (variable_name == "left_pad") {
        return parameters.left_pad;
    } else if (variable_name == "right_pad") {
        return parameters.right_pad;
    } else if (variable_name == "p1") {
        return parameters.coef[0][0];
    } else if (variable_name == "kon1") {
        return parameters.kon[0];
    } else if (variable_name == "keq1") {
        return parameters.keq[0];
    } else if (variable_name == "QE1") {
        return parameters.QE[0];
    } else if (variable_name == "p2") {
        return parameters.coef[1][0];
    } else if (variable_name == "kon2") {
        return parameters.kon[1];
    } else if (variable_name == "keq2") {
        return parameters.keq[1];
    } else if (variable_name == "QE2") {
        return parameters.QE[1];
    } else {
        throw std::runtime_error("Unknown variable name called in variable_location");
    }
}

void
model(const alglib::real_1d_array &control_parameters, alglib::real_1d_array &residuals, int row)
{
    //The 3 Concentration Arrays.
    double** species_in  = parameters.species_in;
    double** species_out = parameters.species_out;
    double* r = parameters.r;

    double model_profile[parameters.row_count * parameters.window_size];

    double left_pad = parameters.left_pad, right_pad = parameters.right_pad;
    /*
    double p1 = parameters.p1, kon1 = parameters.kon1 * parameters.dt, koff1 = parameters.koff1 * parameters.dt, keq1 = parameters.keq1, QE1 = parameters.QE1;
    double p2 = parameters.p2, kon2 = parameters.kon2 * parameters.dt, koff2 = parameters.koff2 * parameters.dt, keq2 = parameters.keq2, QE2 = parameters.QE2;
    double rate_ratio = kon2 / kon1;
    if (parameters.solve_for_kon1) {
        kon2 = rate_ratio * kon1;
        parameters.kon2 = kon2 / parameters.dt;
    }

    koff1 = kon1 / keq1;
    parameters.koff1 = koff1 / parameters.dt;
    koff2 = kon2 / keq2;
    parameters.koff2 = koff2 / parameters.dt;
    */
    
    Eigen::MatrixXd* inverted_diffusion_matrix  = new Eigen::MatrixXd[parameters.number_of_species];
    Eigen::MatrixXd* E                          = new Eigen::MatrixXd[parameters.number_of_species];
    Eigen::MatrixXd* solution                   = new Eigen::MatrixXd[parameters.number_of_species];
    for (int specie = 0; specie < parameters.number_of_species; specie++) {
        E[specie]           = Eigen::MatrixXd::Zero(parameters.X, 1);
        solution[specie]    = Eigen::MatrixXd::Zero(parameters.X, 1);
        inverted_diffusion_matrix[specie] = Eigen::MatrixXd(parameters.X, parameters.X);
        for(int j = 0; j < parameters.X; j++) {
            for (int i = 0; i < parameters.X; i++) {
                inverted_diffusion_matrix[specie](j,i) = (*parameters.inverted_diffusion_matrix)[specie](j,i);
            }
        }
    }
    //double reaction_rate_1 = 0.0d, reaction_rate_2 = 0.0d;
    double scatter = 1.0d;
    double reaction_rate[parameters.number_of_reactions];
    double specie_reaction_rate[parameters.number_of_species], available[parameters.number_of_species];
    for (int specie = 0; specie < parameters.number_of_species; specie++) {
        for (int x = 0; x < parameters.X; x++) {
            solution[specie](x, 0) = species_in[specie][x + row * parameters.X]; // Initialize inlets
        }
    }

    for (int z = 0; z < parameters.Z; z++) {
        for (int i = 0; i < parameters.X; i++) { 
            for (int specie = 0; specie < parameters.number_of_species; specie++) {
                if (i == 0) { // left edge
                    available[specie] =                                            (1.00d - r[specie]) * solution[specie](i, 0)          + r[specie] * solution[specie](i + 1, 0);
                } else if (i == parameters.X - 1) { // right edge
                    available[specie] = r[specie] * solution[specie](i - 1, 0)   + (1.00d - r[specie]) * solution[specie](i, 0);
                } else { // mid points
                    available[specie] = r[specie] * solution[specie](i - 1, 0)   + (2.00d - 2.00d * r[specie]) * solution[specie](i, 0)  + r[specie] * solution[specie](i + 1, 0);
                }
                available[specie] = max(available[specie], 0.0d);
            }
            if (parameters.disable_reactions) {
                for (int reaction = 0; reaction < parameters.number_of_reactions; reaction++) {
                    reaction_rate[reaction] = 0.0d;
                }
                //reaction_rate_1 = 0.0d;
                //reaction_rate_2 = 0.0d;
            } else { // TODO Enable chemical specie / reaction index lookup to remove hard coded values
                //reaction_rate_1 = kon1 * solution[0](i, 0) * solution[1](i, 0) - koff1 * solution[2](i, 0);
                //reaction_rate_2 = kon2 * solution[0](i, 0) * solution[2](i, 0) - koff2 * solution[3](i, 0);
                // specie reaction rates
                
                for (int reaction = 0; reaction < parameters.number_of_reactions; reaction++) {
                    reaction_rate[reaction] = parameters.kon[reaction];
                    double reverse = parameters.kon[reaction] / parameters.keq[reaction];
                    for (int specie = 0; specie < parameters.number_of_species; specie++) {
                        if (parameters.coef[reaction][specie] < 0.0d) {
                            reaction_rate[reaction] *= solution[specie](i, 0); // Forward Reaction
                        } else if (parameters.coef[reaction][specie] > 0.0d) {
                            reverse *= solution[specie](i, 0); // Reverse Reaction
                        }
                    }
                    reaction_rate[reaction] -= reverse;
                }
                // limiting reagents
                for (int specie = 0; specie < parameters.number_of_species; specie++) {
                    double specie_rate = 0.0d;
                    for (int reaction = 0; reaction < parameters.number_of_reactions; reaction++) {
                        specie_rate += parameters.coef[reaction][specie] * reaction_rate[reaction];
                    }
                    while(available[specie] + specie_rate < 0) { // check if NS is limiting reagent.
                        available[specie] = max(available[specie], 0.0d);
                        double total_positive_magnitude = 0.00d;
                        for (int reaction = 0; reaction < parameters.number_of_reactions; reaction++) {
                            if (reaction_rate[reaction] * parameters.coef[reaction][specie] > 0.0d) {
                                total_positive_magnitude += reaction_rate[reaction] * parameters.coef[reaction][specie];
                            }
                        }
                        if (total_positive_magnitude == 0.00d) {
                            for (int reaction = 0; reaction < parameters.number_of_reactions; reaction++) {
                                if (reaction_rate[reaction] * parameters.coef[reaction][specie] < 0.0d) {
                                    reaction_rate[reaction] = 0.0d;
                                }
                            } 
                        } else {
                            for (int reaction = 0; reaction < parameters.number_of_reactions; reaction++) {
                                if (reaction_rate[reaction] * parameters.coef[reaction][specie] < 0) {
                                    if (abs(reaction_rate[reaction] * parameters.coef[reaction][specie]) > total_positive_magnitude) {
                                        reaction_rate[reaction] = std::copysign(total_positive_magnitude / parameters.coef[reaction][specie], reaction_rate[reaction]);
                                    } else {
                                        reaction_rate[reaction] = max(reaction_rate[reaction] * parameters.coef[reaction][specie] * 0.99, reaction_rate[reaction] * parameters.coef[reaction][specie]) / parameters.coef[reaction][specie];
                                    }
                                }  
                            }
                        }
                        specie_rate = 0.0d;
                        for (int reaction = 0; reaction < parameters.number_of_reactions; reaction++) {
                            specie_rate += parameters.coef[reaction][specie] * reaction_rate[reaction];
                            if (std::isnan(specie_rate)) {
                                throw std::runtime_error("On row:" + std::to_string(row) + " z:" + std::to_string(z) + " specie:" + std::to_string(specie) + " i:" + std::to_string(i) + " reaction rate is nan");
                            }
                            if (std::isinf(specie_rate)) {
                                throw std::runtime_error("On row:" + std::to_string(row) + " z:" + std::to_string(z) + " specie:" + std::to_string(specie) + " i:" + std::to_string(i) + " reaction rate is inf");
                            }
                        }
                    }
                    specie_reaction_rate[specie] = specie_rate;
                }
                
                /*
                if(available[1] - reaction_rate_1 / p1 < 0) { // check if NS is limiting reagent.
                    reaction_rate_1 = available[1] * p1;
                }
                if(available[3] + reaction_rate_2 < 0) { // check if NDD is limiting reagent
                    reaction_rate_2 = - available[3];
                }
                while(available[0] - reaction_rate_1 - reaction_rate_2 < 0) { // check if D is limiting reagent
                    reaction_rate_1 = min(reaction_rate_1 * 0.99, reaction_rate_1);
                    reaction_rate_2 = min(reaction_rate_2 * 0.99, reaction_rate_2);
                    if (available[0] == 0.0d && reaction_rate_1 > 0.0d && reaction_rate_2 > 0.0d) {
                        reaction_rate_1 = 0.0d;
                        reaction_rate_2 = 0.0d;
                    } else if (available[0] == 0.0d) {
                        double lesser_magnitude = min(abs(reaction_rate_1),abs(reaction_rate_2));
                        reaction_rate_1 = std::copysign(lesser_magnitude, reaction_rate_1);
                        reaction_rate_2 = std::copysign(lesser_magnitude, reaction_rate_2);
                    }
                }
                while(available[2] + reaction_rate_1 - reaction_rate_2 / p2 < 0) { // check if ND is limiting reagent
                    reaction_rate_1 = max(reaction_rate_1 * 0.99, reaction_rate_1);
                    reaction_rate_2 = min(reaction_rate_2 * 0.99, reaction_rate_2);
                }
            */
            }
            /*
            specie_reaction_rate[0] = -reaction_rate_1 - reaction_rate_2;
            specie_reaction_rate[1] = -reaction_rate_1 / p1;
            specie_reaction_rate[2] = reaction_rate_1 - reaction_rate_2 / p2;
            specie_reaction_rate[3] = reaction_rate_2;
            */
            for (int specie = 0; specie < parameters.number_of_species; specie++) {
                E[specie](i, 0) = available[specie] + specie_reaction_rate[specie];
            }
        }
        for (int specie = 0; specie < parameters.number_of_species; specie++) {
            solution[specie] = inverted_diffusion_matrix[specie] * E[specie];
            for (int i = 0; i < parameters.X; i++) {
                if (std::isnan(solution[specie](i, 0))) {
                    throw std::runtime_error("On row:" + std::to_string(row) + " z:" + std::to_string(z) + " specie:" + std::to_string(specie) + " i:" + std::to_string(i) + " solution is nan");
                }
                if (std::isinf(solution[specie](i, 0))) {
                    throw std::runtime_error("On row:" + std::to_string(row) + " z:" + std::to_string(z) + " specie:" + std::to_string(specie) + " i:" + std::to_string(i) + " solution is inf");
                }
                solution[specie](i, 0) = max(solution[specie](i, 0), 0.00d);

                if ( z == (parameters.Z - 1)) {
                    species_out[specie][i + row * parameters.X] = solution[specie](i, 0);
                }
            }             
        }
    }

    double split = 99.5d;
    int bottom_point, top_point;
    parameters.num_profile_width = (double)parameters.window_size - left_pad - right_pad;    
    for (int i = 0; i < parameters.window_size; i++) {
        if (i < left_pad) {
            model_profile[row * parameters.window_size + i] = 0;
        } else if (i < parameters.window_size - right_pad) {
            split = (double)(i - left_pad) * ((double)parameters.X - 1.0d) / (double)(parameters.num_profile_width - 1.0d);
            bottom_point = static_cast<int>(floor(split));
            top_point = static_cast<int>(ceil(split));
            if (top_point == bottom_point || bottom_point == (parameters.X - 1)) {
                scatter = scattering_correction(species_out[1][bottom_point + row * parameters.X], species_out[2][bottom_point + row * parameters.X], parameters.coef[0][0]);
                for (int specie = 0; specie < parameters.number_of_species; specie++) {
                    model_profile[row * parameters.window_size + i] += species_out[specie][bottom_point + row * parameters.X] * parameters.QE[specie] * scatter;
                }
            } else {
                scatter = scattering_correction(species_out[1][bottom_point + row * parameters.X] * ((double)top_point - split) + species_out[1][top_point + row * parameters.X] * (split - (double)bottom_point), species_out[2][bottom_point + row * parameters.X] * ((double)top_point - split) + species_out[2][top_point + row * parameters.X] * (split - (double)bottom_point), parameters.coef[0][0]);
                for (int specie = 0; specie < parameters.number_of_species; specie++) {
                    model_profile[row * parameters.window_size + i] += ((species_out[specie][bottom_point + row * parameters.X] * parameters.QE[specie] ) * ((double)top_point - split) + (species_out[specie][top_point + row * parameters.X] * parameters.QE[specie]) * (split - (double)bottom_point)) * scatter;
                }
            }
        } else {
            model_profile[row * parameters.window_size + i] = parameters.dye_conc[row];
        }
        residuals[i + row * parameters.window_size] = pow(model_profile[row * parameters.window_size + i] / parameters.dye_conc[row] - parameters.experimental_profile[row * parameters.window_size + i], 2.0d);
        parameters.numeric_model_profile[row * parameters.window_size + i] = model_profile[row * parameters.window_size + i];
    }

    delete[] E;
    delete[] solution;
    delete[] inverted_diffusion_matrix;
}



void
alglib_solver(const alglib::real_1d_array &control_parameters, alglib::real_1d_array &residuals, void *ptr)
{
    std::cout << parameters.iterations;
    if (!parameters.run_solver) {
        return;
    } else {
        for (int index = 0; index < parameters.number_of_variables; index++) { // This sorts out out of order parameters
            variable_location(parameters.solve_for[index]) = control_parameters[index];
        }
        thread* threads = new thread[parameters.row_count];
        for (int row = 0; row < parameters.row_count; row++) {
            threads[row] = thread(model, control_parameters, std::ref(residuals), row);
        }
        for (int row = 0; row < parameters.row_count; row++) {
            threads[row].join();
        }
        delete[] threads;
    }
    for (int i = 0; i < std::to_string(parameters.iterations).length(); i++) {
        std::cout << '\b' << ' ' << '\b';
    }
    parameters.iterations = parameters.iterations + 1;
}

int
main()
{
    std::cout << "Starting main" << endl;

    std::cout << "Openning Database navier.db: ";
    sqlite3 *db; // Database connection handle
    int rc = sqlite3_open("../../navier.db", &db);
    char* zErrMsg = 0; // For error messages
    if (rc != SQLITE_OK) {
        fprintf(stderr, "Cannot open database: %s\n", sqlite3_errmsg(db));
    } else {
        fprintf(stdout, "Database opened successfully\n");
    }
    std::cout << "Retrieving model control parameters:";
    read_model_parameters_from_db(db);
    parameters.file_name = "../" + parameters.experiment_name;
    std::cout << " Done" << endl;

    std::cout << "Retrieving line count:";
    lines_from_profile_text(db);
    std::cout << " Done" << endl;

    std::cout << "Reading experimental parameters from database: ";
    parameters.wt_percent                = new double[parameters.row_count]();   //= 0.1d; // wt%
    parameters.dye_conc_mgml             = new double[parameters.row_count]();   //= 0.00336d; // mg/ml FITC 
    parameters.dye_conc                  = new double[parameters.row_count]();   //= dye_conc_mgml * 1000.0d / 332.326d * 6.022e+23;     // molecules FITC / m3
    parameters.beginning_of_channel      = new int[parameters.row_count]();
    parameters.end_of_channel            = new int[parameters.row_count]();
    read_exp_parameters_from_db(db);
    std::cout << "Done" << endl;

    std::cout << parameters.number_of_variables << " Variables:(";
    for (int item = 0; item < parameters.number_of_variables; item++) {
        std::cout << parameters.solve_for[item];
        if (item != parameters.number_of_variables - 1) {std::cout << " ";} else {std::cout << ")" << endl;}
    }
    std::cout <<  parameters.number_of_species << " Species:(";
    for (int item = 0; item < parameters.number_of_species; item++) {
        std::cout << parameters.SPECIE[item];
        if (item != parameters.number_of_species - 1) {std::cout << " ";} else {std::cout << ")"<< endl;}
    }
    std::cout << parameters.number_of_reactions << " Reactions:(";
    for (int item = 0; item < parameters.number_of_reactions; item++) {
        std::cout << parameters.reaction[item];
        if (item != parameters.number_of_reactions - 1) {std::cout << " ";} else {std::cout << ")"<< endl;}
    }
    parameters.output_file_name = parameters.file_name + " " + parameters.normalization_method + " scatter_" + parameters.scatter_correction_type + " output.txt";

    parameters.CONC = new double**[parameters.number_of_entrances];
    for (int entrance = 0; entrance < parameters.number_of_entrances; entrance++) {
        parameters.CONC[entrance] = new double*[parameters.number_of_species];
        for (int specie = 0; specie < parameters.number_of_species; specie++) {
            parameters.CONC[entrance][specie] = new double[parameters.row_count]();
        }
    }
    
    std::cout << "Reading specie and reaction values from database: "; // TODO introduce species lookup and unit correction
    parameters.QE                      = new double[parameters.number_of_species];
    parameters.diffusion_rate          = new double[parameters.number_of_species];
    parameters.diameter                = new double[parameters.number_of_species];
    parameters.species_unit_conversion = new double[parameters.number_of_species]();
    for (int specie = 0; specie < parameters.number_of_species; specie++) {
        parameters.specie_type.push_back("particle"); // Establishes size of specie_type vector
        if(specie == 1) {
            parameters.species_unit_conversion[specie] = 1.00d / 100.0d / 1.05d / ( 4.0d / 3.0d * M_PI * pow(parameters.diameter_old / 2.0d, 3.0d)) / 1000 / (6.022e+23) * 1.0e+6; // converts units between umolar beads and wt% beads
        } else {
            parameters.species_unit_conversion[specie] = 332.326d / 1.0e+6;
        }
    }
    parameters.kon                     = new double[parameters.number_of_reactions];
    parameters.keq                     = new double[parameters.number_of_reactions];
    parameters.coef                    = new double*[parameters.number_of_reactions];
    parameters.exp                     = new double*[parameters.number_of_reactions];
    for (int react = 0; react < parameters.number_of_reactions; react++) {
        parameters.react_ks_text.push_back("foo"); //Establishes size of reaction vectors
        parameters.react_exp_text.push_back("foo"); 
        parameters.react_coef_text.push_back("foo");
        parameters.react_specie_text.push_back("foo");
        parameters.coef[react]         = new double[parameters.number_of_species];
        parameters.exp[react]          = new double[parameters.number_of_species];
    }
    read_specie_and_reaction_values_from_db(db);
    double total_flowrate = 0.0d;
    for (int entrance = 0; entrance < parameters.number_of_entrances; entrance++) {
        total_flowrate += parameters.ENTRANCE_FLOWRATE[entrance];
    }
    double restime          = parameters.W * parameters.H * parameters.L / (total_flowrate); // seconds
    parameters.dt           = restime / parameters.Z; // seconds
    parameters.r            = new double[parameters.number_of_species]();
    Eigen::MatrixXd* inverted_diffusion_matrix  = new Eigen::MatrixXd[parameters.number_of_species];
    parameters.inverted_diffusion_matrix  = &inverted_diffusion_matrix;
    for (int specie = 0; specie < parameters.number_of_species; specie++) {
        inverted_diffusion_matrix[specie] = Eigen::MatrixXd::Zero(parameters.X, parameters.X);
        if (parameters.specie_type[specie] == "particle" && parameters.diameter[specie] != 0.0d) {
            parameters.diffusion_rate[specie] = 1.380649e-23d * parameters.temperature / (3.0d * M_PI * parameters.visc * parameters.diameter[specie]);
        }
        parameters.r[specie] = parameters.diffusion_rate[specie] * parameters.dt / (parameters.W * parameters.W) * parameters.X * parameters.X; // r = dT/(dX)^2, dT = dif * dt / W^2, 1/dX = X, T = Dt/l^2
        for (int i = 0; i < parameters.X; i++) {
            for (int j = 0; j < parameters.X; j++) {
                if (j == i - 1 || j == i + 1) {
                    inverted_diffusion_matrix[specie](i, j) = -parameters.r[specie];
                } else if ((i != 0 && i != parameters.X - 1) && (i == j)) {
                    inverted_diffusion_matrix[specie](i, j) = 2.0d + 2.0d * parameters.r[specie];
                } else if ((i == 0 || i == parameters.X - 1) && (i == j)) {
                    inverted_diffusion_matrix[specie](i, j) = 1.0d + parameters.r[specie];
                } else {
                    inverted_diffusion_matrix[specie](i, j) = 0.0d;
                }
            }
        }
        inverted_diffusion_matrix[specie] = inverted_diffusion_matrix[specie].inverse();
    }
    std::cout << "Done" << endl;

    for (int react = 0; react < parameters.number_of_reactions; react++) {
        parameters.kon[react] *= parameters.dt;
        std::cout << "Reaction:" << parameters.reaction[react] << " ";
        for (int specie = 0; specie < parameters.number_of_species; specie++) {
            std::cout << parameters.coef[react][specie] << " " << parameters.SPECIE[specie];
            if (specie != parameters.number_of_species - 1) {std::cout << ", ";}
        }
        std::cout << endl;
    }

    std::cout << "Reading experimental profiles from database: ";
    parameters.experimental_profile    = new double[parameters.row_count * parameters.window_size]{};
    parameters.numeric_model_profile   = new double[parameters.row_count * parameters.window_size]{};
    read_raw_profiles_from_db(db);
    // The inlet and outlet Concentration Arrays.
    parameters.species_in    = new double*[parameters.number_of_species];
    parameters.species_out   = new double*[parameters.number_of_species];
    for (int specie = 0; specie < parameters.number_of_species; specie++) {
        parameters.species_in[specie]    = new double[parameters.X * parameters.row_count]{};
        parameters.species_out[specie]   = new double[parameters.X * parameters.row_count]{};
    } 

    for (int row = 0; row < parameters.row_count; row++) {
        //initialize the concentration inlet arrays. 
        parameters.dye_conc_mgml[row] = max(parameters.CONC[0][0][row],parameters.CONC[1][0][row]);
        parameters.dye_conc[row] = parameters.dye_conc_mgml[row] / parameters.species_unit_conversion[0];
        double FLOWRATE_ACCOUNTED_FOR = 0.0d;
        for (int entrance = static_cast<int>(parameters.X * FLOWRATE_ACCOUNTED_FOR / total_flowrate); entrance < parameters.number_of_entrances; entrance++) {
            for (int x = 0; x < parameters.X; x++) {
                int i = x + parameters.X * row;
                for (int specie = 0; specie < parameters.number_of_species; specie++) {
                    if (x < (static_cast<int>(parameters.X * parameters.ENTRANCE_FLOWRATE[entrance] / total_flowrate))) {
                        parameters.species_in[specie][i] = parameters.CONC[entrance][specie][row] / parameters.species_unit_conversion[specie];
                    }
                }
            }
            FLOWRATE_ACCOUNTED_FOR += parameters.ENTRANCE_FLOWRATE[entrance];
        }
    }
    std::cout << "Done" << endl;

    std::cout << "Reading alglib inputs from database: ";
    parameters.scale                    = new double[parameters.number_of_variables];
    parameters.initial_values_alglib    = new double[parameters.number_of_variables];
    parameters.low_bound                = new double[parameters.number_of_variables];
    parameters.up_bound                 = new double[parameters.number_of_variables];
    read_alglib_values_from_db(db);
    std::cout << "Done" << endl;

    std::cout << "Normalizing " << parameters.row_count << " profiles with " << parameters.normalization_method << ": ";
    normalize_profile();
    std::cout << "Done" << endl;
    
    if (parameters.save_normalized_profiles) {
        std::cout << "Writing " << parameters.row_count << " normalized profiles to database: ";
        write_normalized_values_to_db(db);
        std::cout << "Done" << endl;
    }

    std::cout << "Priming " << parameters.output_file_name << ": ";
    prime_excel_output(parameters.output_file_name);
    std::cout << "Done" << endl;

    try
    {
        /*
        double epsx = 1e-11;
        finishes if |v|<=EpsX is fulfilled 
        |.| means Euclidian norm
        v - scaled step vector
        v[i]=dx[i]/s[i]
        dx - step vector
        dx=X(k+1)-X(k)
        s - scaling coefficients set by MinLMSetScale()
        Recommended values: 1E-9 ... 1E-12.
        */
        double DiffStep = 0.0001;
        
        alglib::real_1d_array control_parameters;
        control_parameters.setcontent(parameters.number_of_variables, parameters.initial_values_alglib);
        alglib::real_1d_array s;
        s.setcontent(parameters.number_of_variables, parameters.scale);
        alglib::real_1d_array bndl;
        bndl.setcontent(parameters.number_of_variables, parameters.low_bound);
        alglib::real_1d_array bndu;
        bndu.setcontent(parameters.number_of_variables, parameters.up_bound);
        alglib::minlmstate state;
        alglib::minlmreport rep;
        alglib::minlmcreatev(parameters.number_of_variables, parameters.window_size * parameters.row_count, control_parameters, DiffStep, state);
        alglib::minlmsetbc(state, bndl, bndu);
        alglib::minlmsetcond(state, parameters.convergence_epsx, parameters.max_iterations);
        alglib::minlmsetscale(state, s);
        alglib::minlmsetnonmonotonicsteps(state, 2);
        if (parameters.run_solver) {
            std::cout << "minlmoptimize: " << parameters.experiment_name << " on iteration: ";
            alglib::minlmoptimize(state, alglib_solver);   // Optimize
            std::cout << parameters.iterations << ". Done" << endl;
            alglib::minlmresults(state, control_parameters, rep);
            printf("%s\n", control_parameters.tostring(4).c_str());
        };

        for (int index = 0; index < parameters.number_of_variables; index++) { // This sorts out of order parameters
            parameters.initial_values_alglib[index] = control_parameters[index];
            variable_location(parameters.solve_for[index]) = control_parameters[index];
        }
        if (parameters.solve_for_keq1) {parameters.koff1 = parameters.kon1 / parameters.keq1;}
        if (parameters.solve_for_keq2) {parameters.koff2 = parameters.kon2 / parameters.keq2;}
 
        save_excel_output();

        if (parameters.save_model_profiles) {
            std::cout << "Writing model profiles to database: ";
            write_model_profile_to_db(db);
            std::cout << "Done" << endl;
        }

        string accept_alglib_values;
        cout << "Save parameter solutions to db as initial values? (y/n) ";
        cin >> accept_alglib_values;
        if (accept_alglib_values == "y") {
            std::cout << "Writing alglib_values to database: ";
            write_alglib_values_to_db(db);
            std::cout << "Done" << endl;
        }
        
    }
    catch(alglib::ap_error alglib_exception)
    {
        printf("ALGLIB exception with message '%s'\n", alglib_exception.msg.c_str());
        return 1;
    }

    sqlite3_close(db);
    std::cout << "Database closed." << endl;

    std::cout << "Program finished " << endl;
    return 0;
}