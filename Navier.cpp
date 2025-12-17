#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include <eigen-3.4.0/Eigen/Dense>
#include <optimization.h>
#include <stdafx.h>
#include <solvers.h>
#include <sqlite3.h>
// D: Dye
// NS: Nanoparticle Site
// ND: Nanoparticle Bound Dye

using namespace std;
// Text FILE HANDLERS
void prime_excel_output(string file_name);
void save_excel_output(double* D, double* ND, double* NS, double* species_4);

// DB HANDLERS
void lines_from_profile_text(sqlite3* db);
void delete_values_from_db(sqlite3* db, string table, string where_conditions);
int exp_parameters_db_callback(void *data, int count, char **argv, char **columnNames);
void read_exp_parameters_from_db(sqlite3* db);
int raw_profiles_db_callback(void *data, int count, char **argv, char **columnNames);
void read_raw_profiles_from_db(sqlite3* db);
int alglib_input_db_callback(void *data, int count, char **argv, char **columnNames);
void read_alglib_values_from_db(sqlite3* db);
void write_normalized_values_to_db(sqlite3* db);
void write_model_profile_to_db(sqlite3* db, double* D, double* ND, double* NS, double* species_4);

// MODEL
void normalize_profile();
void solver(double* control_parameters, alglib::real_1d_array &residuals);
void alglib_solver(const alglib::real_1d_array &control_parameters, alglib::real_1d_array &residuals, void *ptr);
double scattering_correction(double NS, double ND, double p);

struct parameters_struct {
    // Device Dimenssions
    int Z, X, M;
    const double W = 5e-4, H = 4e-5, L = 0.025;  //meters: 500 um, 40 um, 2.5 cm 

    // Operating conditions/settings
    string normalization_method, scatter_correction_type;
    bool run_solver = true, disable_reactions = false;// = false;
    string experiment_name, file_name, output_file_name;// = "../" + experiment_name;

    double ENTRANCE_1_FLOWRATE;
    double ENTRANCE_2_FLOWRATE;
    string ENTRANCE_1_SPECIES_1_NAME;
    double* ENTRANCE_1_SPECIES_1_CONC;
    string ENTRANCE_2_SPECIES_1_NAME;
    double* ENTRANCE_2_SPECIES_1_CONC;
    string ENTRANCE_1_SPECIES_2_NAME;
    double* ENTRANCE_1_SPECIES_2_CONC;
    string ENTRANCE_2_SPECIES_2_NAME;
    double* ENTRANCE_2_SPECIES_2_CONC;

    double flow;                            // m3/s: 2*5 ulmin
    double diameter;                        // meters
    double* wt_percent;                     //= 0.1d; // wt%
    double* dye_conc_mgml;                  //= 0.00336d; // mg/ml FITC
    double* bead_conc;                      //= wt_percent / 100.0d / 1.05d / ( 4.0d / 3.0d * M_PI * pow(diameter / 2.0d, 3.0d)); // beads/m3 
    double* dye_conc;                       //= dye_conc_mgml * 1000.0d / 332.326d * 6.022e+23;     // molecules FITC / m3
    int* beginning_of_channel;
    int* end_of_channel;
    double experimental_profile[31][180];
    double numeric_model_profile[31][180];

    double* D_in;
    double* NS_in;
    double* ND_in;
    double* species_4_in;
    double* D;
    double* NS;
    double* ND;
    double* species_4;
    double* D_out;
    double* NS_out;
    double* ND_out;
    double* species_4_out;

    int row_count = 0, window_size = 180, exp_left_padding = 20, exp_right_padding = 20;
    double num_profile_width;
    
    int low_ref_start, low_ref_end, high_ref_start, high_ref_end, iterations = 0;

    // Physical Constants
    double visc; // Dynamic viscosity of water at 20C in Pa.s
    double difusion_dye; // m2/s From 4.9 × 10−6 cm2 s−1 The diffusion coefficient of fluorescein in water at 21.5°C, as calculated from the Wilke-Chang correlation
    double difusion_beads; // m2/s From kB*T/(3*pi*visc*d) kB=1.380649×10−23 J⋅K−1

    double restime; // seconds
    double dt; // seconds
    double T_dye;
    double T_beads;
    double dT_dye; // T = Dt/l^2
    double dT_beads;
    double dx; // m / subdivision
    double dX; // X = x/l
    double r_dye; // r = dT/(dX)^2
    double r_beads; // r = dT/(dX)^2

    int number_of_variables;
    double left_pad, right_pad, p1, kon1, koff1, keq1, QE1, p2, kon2, koff2, keq2, QE2;

    string PARAMETERS_TO_SOLVE_FOR;
    vector<string> solve_for;
    bool solve_for_left_pad = false, solve_for_right_pad = false;
    bool solve_for_p1 = false, solve_for_kon1 = false, solve_for_koff1 = false, solve_for_keq1 = false, solve_for_QE1 = false;
    bool solve_for_p2 = false, solve_for_kon2 = false, solve_for_koff2 = false, solve_for_keq2 = false, solve_for_QE2 = false;

    double* initial_values_alglib;
    double* low_bound;
    double* up_bound; 
    double* scale;
} parameters;


struct matrix_struct {
    Eigen::Matrix<double, -1, -1>* m;
    Eigen::Matrix<double, -1, -1>* m_beads;
    Eigen::Matrix<double, -1, -1>* E_D;
    Eigen::Matrix<double, -1, -1>* E_NS;
    Eigen::Matrix<double, -1, -1>* E_ND;
    Eigen::Matrix<double, -1, -1>* E_species_4;
    Eigen::Matrix<double, -1, -1>* solution_D;
    Eigen::Matrix<double, -1, -1>* solution_NS;
    Eigen::Matrix<double, -1, -1>* solution_ND;
    Eigen::Matrix<double, -1, -1>* solution_species_4;
} matrixes;

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
save_excel_output(double* D, double* ND, double* NS, double* species_4)
{
    double dye_bead_ratio = 2000.0d;
    double total_dye = 0.0d;
    int z = parameters.Z - 1;
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
            num_derivative[i] = (parameters.numeric_model_profile[row][i+1] / parameters.dye_conc[row] - parameters.numeric_model_profile[row][i-1] / parameters.dye_conc[row]) / (2 * parameters.W * 1.0e+6 / parameters.window_size);
            exp_derivative[i] = (parameters.experimental_profile[row][i+1] - parameters.experimental_profile[row][i-1]) / (2 * parameters.W * 1.0e+6 / parameters.window_size);
        }
        for (int j = 0; j < 10; j++) {
            fout << z * parameters.dt << "  " << parameters.p1 << "  " << parameters.kon1 << "  " << parameters.koff1 << "  " << parameters.dye_conc[row] << "   " <<  parameters.bead_conc[row] << "  " <<  parameters.QE1 << "  ";
            switch(j) {
                case 0:
                    fout << "Free_Dye  ";
                    for (int x = 0; x < parameters.X; x++) {
                        fout << D[x + row * parameters.X] / parameters.dye_conc[row] << "    ";
                    }
                    break;
                case 1:
                    fout << "Bound_Dye  ";
                    for (int x = 0; x < parameters.X; x++) {
                        fout << (ND[x + row * parameters.X] + species_4[x + row * parameters.X]) / parameters.dye_conc[row] << "    ";
                    }
                    break;
                case 2:
                    fout << "Total_Dye  ";
                    for (int x = 0; x < parameters.X; x++) {
                        fout << (D[x + row * parameters.X] + ND[x + row * parameters.X] + species_4[x + row * parameters.X]) / parameters.dye_conc[row] << "    ";
                    }
                    break;
                case 3:
                    fout << "Unbound_Beads  ";
                    for (int x = 0; x < parameters.X; x++) {
                        fout << NS[x + row * parameters.X] * 100.0d * 1.05d * ( 4.0d / 3.0d * M_PI * pow(parameters.diameter / 2.0d, 3.0d)) << "    ";

                    }
                    break;
                case 4:
                    fout << "Bound_Beads  ";
                    for (int x = 0; x < parameters.X; x++) {
                        fout << (ND[x + row * parameters.X] / parameters.p1 + species_4[x + row * parameters.X] / parameters.p1 / parameters.p2)* 100.0d * 1.05d * ( 4.0d / 3.0d * M_PI * pow(parameters.diameter / 2.0d, 3.0d)) << "    ";
                    }
                    break;   
                case 5:
                    fout << "Total_Beads  ";
                    for (int x = 0; x < parameters.X; x++) {
                        fout << (NS[x + row * parameters.X] + ND[x + row * parameters.X] / parameters.p1 + species_4[x + row * parameters.X] / parameters.p1 / parameters.p2) * 100.0d * 1.05d * ( 4.0d / 3.0d * M_PI * pow(parameters.diameter / 2.0d, 3.0d)) << "    ";
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
                        fout << parameters.experimental_profile[row][x] << "    ";
                    }
                    break;
                case 9:
                    fout << "Total_Dye_rescale  ";
                    for (int x = 0; x < parameters.window_size; x++) {
                        fout << parameters.numeric_model_profile[row][x] / parameters.dye_conc[row] << "    ";
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
        if (criterion == "WT_PERCENT") {
            if (argv[i] != NULL) {
                parameters.row_count++;
            } else if (argv[i] == NULL) {
                throw 101;
            }
        }
    };
    return 0;
}


void
lines_from_profile_text(sqlite3* db)
{
    char* errMsg = 0;
    string sqltext = "SELECT WT_PERCENT FROM 'raw_profile' WHERE NAME='" + parameters.experiment_name + "';";
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
    string where_conditions_used = "1=2";
    if(!where_conditions.empty()) {
        where_conditions_used = where_conditions;
    }
    sqltext = sqltext + "DELETE FROM " + table + " WHERE " + where_conditions_used + ";";
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
        if (!value.empty()) {
            parameters.X = stoi(value);   
        } else if (value.empty()) {
            parameters.X = 500;
        }
    } else if (criterion == "length/time resolution (Z)") {
        if (!value.empty()) {
            parameters.Z = stoi(value);
        } else if (value.empty()) {
            parameters.Z = 100;
        }
    } else if (criterion == "experiment_name") {
        if (!value.empty()) {
            parameters.experiment_name = value;
        } else if (value.empty()) {
            parameters.experiment_name = "24_3";
        }
    } else if (criterion == "exp_left_padding") {
        if (!value.empty()) {
            parameters.exp_left_padding = stoi(value);
        } else if (value.empty()) {
            parameters.exp_left_padding = 20;
        }
    } else if (criterion == "exp_right_padding") {
        if (!value.empty()) {
            parameters.exp_right_padding = stoi(value);
        } else if (value.empty()) {
            parameters.exp_right_padding = 20;
        }
    } else if (criterion == "disable_reactions") {
        if (!value.empty()) {
            if(value == "true") {
                parameters.disable_reactions = true;
            } else {
                parameters.disable_reactions = false;
            }
        } else if (value.empty()) {
            parameters.disable_reactions = false;
        }
    } else if (criterion == "run_solver") {
        if (!value.empty()) {
            if(value == "true") {
                parameters.run_solver = true;
            } else {
                parameters.run_solver = false;
            }
        } else if (value.empty()) {
            parameters.run_solver = true;
        }
    } else if (criterion == "scatter_correction_type") {
        if (!value.empty()) {
            parameters.scatter_correction_type = value;
        } else if (value.empty()) {
            parameters.scatter_correction_type = "none";
        }
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
                    parameters.dye_conc[row] = parameters.dye_conc_mgml[row] * 1000.0d / 332.326d * 6.022e+23;   
                } else if (argv[i] == NULL) {
                    parameters.dye_conc_mgml[row] = 0.0d;
                    parameters.dye_conc[row] = 0.0d;
                }
            }
        } else if (criterion == "BEAD_DIAMETER_NM") {
            if (argv[i] != NULL) {
                parameters.diameter = stod(argv[i]) * 1.0e-9d;
            } else if (argv[i] == NULL) {
                parameters.diameter = 41.0e-9d;
            }
        } else if (criterion == "FLOW_RATE") {
            if (argv[i] != NULL) {
                parameters.flow = stod(argv[i]);
            } else if (argv[i] == NULL) {
                parameters.flow = 0.00000000016667;
            }
        } else if (criterion == "LOW_REF_LEFT") {
            if (argv[i] != NULL) {
                parameters.low_ref_start = stoi(argv[i]);
            } else if (argv[i] == NULL) {
                parameters.low_ref_start = 0;
            }
        } else if (criterion == "LOW_REF_RIGHT") {
            if (argv[i] != NULL) {
                parameters.low_ref_end = stoi(argv[i]);
            } else if (argv[i] == NULL) {
                parameters.low_ref_end = 10;
            }
        } else if (criterion == "HIGH_REF_LEFT") {
            if (argv[i] != NULL) {
                parameters.high_ref_start = stoi(argv[i]);
            } else if (argv[i] == NULL) {
                parameters.high_ref_start = 95;
            }
        } else if (criterion == "HIGH_REF_RIGHT") {
            if (argv[i] != NULL) {
                parameters.high_ref_end = stoi(argv[i]);
            } else if (argv[i] == NULL) {
                parameters.high_ref_end = 105;
            }
        } else if (criterion == "ENTRANCE_1_FLOWRATE") {
            if (argv[i] != NULL) {
                parameters.ENTRANCE_1_FLOWRATE = stod(argv[i]);
            } else if (argv[i] == NULL) {
                parameters.ENTRANCE_1_FLOWRATE = 0.00000000016667 / 2.0d;
            }
        } else if (criterion == "ENTRANCE_2_FLOWRATE") {
            if (argv[i] != NULL) {
                parameters.ENTRANCE_2_FLOWRATE = stod(argv[i]);
            } else if (argv[i] == NULL) {
                parameters.ENTRANCE_2_FLOWRATE = 0.00000000016667 / 2.0d;
            }
        } else if (criterion == "ENTRANCE_1_SPECIES_1_NAME") {
            if (argv[i] != NULL) {
                parameters.ENTRANCE_1_SPECIES_1_NAME = argv[i];
            } else if (argv[i] == NULL) {
            }
        } else if (criterion == "ENTRANCE_1_SPECIES_1_CONC") {
            for (int row = 0; row < parameters.row_count; row++) {
                if (argv[i] != NULL) {
                    parameters.ENTRANCE_1_SPECIES_1_CONC[row] = stod(argv[i]);
                } else if (argv[i] == NULL) {
                    parameters.ENTRANCE_1_SPECIES_1_CONC[row] = 0.0d;
                }
            }
        } else if (criterion == "ENTRANCE_2_SPECIES_1_NAME") {
            if (argv[i] != NULL) {
                parameters.ENTRANCE_2_SPECIES_1_NAME = argv[i];
            } else if (argv[i] == NULL) {
            }
        } else if (criterion == "ENTRANCE_2_SPECIES_1_CONC") {
            for (int row = 0; row < parameters.row_count; row++) {
                if (argv[i] != NULL) {
                    parameters.ENTRANCE_2_SPECIES_1_CONC[row] = stod(argv[i]);
                } else if (argv[i] == NULL) {
                    parameters.ENTRANCE_2_SPECIES_1_CONC[row] = 0.0d;
                }
            }
        } else if (criterion == "ENTRANCE_1_SPECIES_2_NAME") {
            if (argv[i] != NULL) {
                parameters.ENTRANCE_1_SPECIES_2_NAME = argv[i];
            } else if (argv[i] == NULL) {
            }
        } else if (criterion == "ENTRANCE_1_SPECIES_2_CONC") {
            for (int row = 0; row < parameters.row_count; row++) {
                if (argv[i] != NULL) {
                    parameters.ENTRANCE_1_SPECIES_2_CONC[row] = stod(argv[i]);
                } else if (argv[i] == NULL) {
                    parameters.ENTRANCE_1_SPECIES_2_CONC[row] = 0.0d;
                }
            }
        } else if (criterion == "ENTRANCE_2_SPECIES_2_NAME") {
            if (argv[i] != NULL) {
                parameters.ENTRANCE_2_SPECIES_2_NAME = argv[i];
            } else if (argv[i] == NULL) {
            }
        } else if (criterion == "ENTRANCE_2_SPECIES_2_CONC") {
            for (int row = 0; row < parameters.row_count; row++) {
                if (argv[i] != NULL) {
                    parameters.ENTRANCE_2_SPECIES_2_CONC[row] = stod(argv[i]);
                } else if (argv[i] == NULL) {
                    parameters.ENTRANCE_2_SPECIES_2_CONC[row] = 0.0d;
                }
            }
        }
        else if (criterion == "PARAMETERS_TO_SOLVE_FOR") {
            parameters.PARAMETERS_TO_SOLVE_FOR = argv[i];
            string s; // variable to store token obtained from the original string
            stringstream ss(parameters.PARAMETERS_TO_SOLVE_FOR); // constructing stream from the string
            while (getline(ss, s, ' ')) { 
                parameters.solve_for.push_back(s); // store token string in the vector
            }
            parameters.number_of_variables = parameters.solve_for.size();
            std::cout << parameters.number_of_variables << " Variables: (" << parameters.PARAMETERS_TO_SOLVE_FOR << ") ";
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
                } else if (argv[i] == NULL) {

                }
            }
        } else if (criterion == "DEFAULT_NORMALIZATION") {
            if (argv[i] != NULL) {
                parameters.normalization_method = argv[i];
            } else if (argv[i] == NULL) {
                parameters.normalization_method = "none";
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
                parameters.bead_conc[row] = parameters.wt_percent[row] / 100.0d / 1.05d / ( 4.0d / 3.0d * M_PI * pow(parameters.diameter / 2.0d, 3.0d));
            } else if (argv[i] == NULL) {
                throw 101;
            }
        } else if (criterion == "CHANNEL_LEFT_EDGE") {
            if (argv[i] != NULL) {
                parameters.beginning_of_channel[row] = stoi(argv[i]);
            } else if (argv[i] == NULL) {
                parameters.beginning_of_channel[row] = 57;
            }
        } else if (criterion == "CHANNEL_RIGHT_EDGE") {
            if (argv[i] != NULL) {
                parameters.end_of_channel[row] = stoi(argv[i]);
            } else if (argv[i] == NULL) {
                parameters.end_of_channel[row] = 168;
            }
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
                        if (col_index < 0) { // Ignore these data points

                        } else if (col_index < parameters.exp_left_padding) { // Points added in left padding
                            parameters.experimental_profile[row][col_index] = 0.0;
                        } else if (col_index < parameters.end_of_channel[row] + parameters.beginning_of_channel[row] - parameters.exp_left_padding) { // Actual data 
                            parameters.experimental_profile[row][col_index] = stod(intensity[x]);
                            last_profile_point = stod(intensity[x]);
                        } else if (col_index < parameters.end_of_channel[row] + parameters.beginning_of_channel[row] - parameters.exp_left_padding) { // Points added in right padding
                            parameters.experimental_profile[row][col_index] = last_profile_point;
                        }
                    } else if (argv[i] == NULL) {
                        parameters.experimental_profile[row][col_index] = 1;
                    }
                }
            }
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
    string s; // variable to store token obtained from the original string
    stringstream ss(parameters.PARAMETERS_TO_SOLVE_FOR); // constructing stream from the string
    while (getline(ss, s, ' ')) { 
        parameters.solve_for.push_back(s); // store token string in the vector
    }
    for(int i = 0; i < count; i++) {
        criterion = columnNames[i];
        if (criterion == "VARIABLE") {
            if (argv[i] != NULL) {
                current_variable = argv[i];
                varible_to_solve_for = (parameters.solve_for_left_pad == true && current_variable == "left_pad") || (parameters.solve_for_right_pad == true && current_variable == "right_pad" ) || (parameters.solve_for_p1 == true && current_variable == "p1") || (parameters.solve_for_kon1 == true && current_variable == "kon1") || (parameters.solve_for_koff1 == true && current_variable == "koff1") || (parameters.solve_for_keq1 == true && current_variable == "keq1") || (parameters.solve_for_QE1 == true && current_variable == "QE1") || (parameters.solve_for_p2 == true && current_variable == "p2") || (parameters.solve_for_kon2 == true && current_variable == "kon2") || (parameters.solve_for_koff2 == true && current_variable == "koff2") || (parameters.solve_for_keq2 == true && current_variable == "keq2") || (parameters.solve_for_QE2 == true && current_variable == "QE2");
                
                for (int item = 0; item < parameters.number_of_variables; item++) {
                    if (parameters.solve_for[item] == current_variable) {
                        index = item;
                    }
                }
            } else if (argv[i] == NULL) {
                throw 101;
            }
        } else if (criterion == "INITIAL VALUE") {
            if (argv[i] != NULL) {
                if (current_variable      == "left_pad")    {parameters.left_pad = stod(argv[i]);}
                else if(current_variable  == "right_pad")   {parameters.right_pad = stod(argv[i]);}
                else if (current_variable == "p1")          {parameters.p1 = stod(argv[i]);}
                else if (current_variable == "kon1")        {parameters.kon1 = stod(argv[i]);}
                else if (current_variable == "koff1")       {parameters.koff1 = stod(argv[i]);}
                else if (current_variable == "QE1")         {parameters.QE1 = stod(argv[i]);}
                else if (current_variable == "p2")          {parameters.p2 = stod(argv[i]);}
                else if (current_variable == "kon2")        {parameters.kon2 = stod(argv[i]);}
                else if (current_variable == "koff2")       {parameters.koff2 = stod(argv[i]);}
                else if (current_variable == "QE2")         {parameters.QE2 = stod(argv[i]);}
                else if (current_variable == "keq1")        {parameters.keq1 = stod(argv[i]);}
                else if (current_variable == "keq2")        {parameters.keq2 = stod(argv[i]);}
                if (varible_to_solve_for) {
                    parameters.initial_values_alglib[index] = stod(argv[i]); // initial values for alglib solver
                };
            } else if (argv[i] == NULL) {
                throw 102;
            }
        } else if (criterion == "LOWER BOUND") {
            if (argv[i] != NULL) {
                if (varible_to_solve_for) {
                    parameters.low_bound[index] = stod(argv[i]);
                }
            } else if (argv[i] == NULL) {
                throw 102;
            }
        } else if (criterion == "UPPER BOUND") {
            if (argv[i] != NULL) {
                if (varible_to_solve_for) {
                    parameters.up_bound[index] = stod(argv[i]);
                }
            } else if (argv[i] == NULL) {
                throw 102;
            }
        } else if (criterion == "SCALE") {
            if (argv[i] != NULL) {
                if (varible_to_solve_for) {
                    parameters.scale[index] = stod(argv[i]);
                }
            } else if (argv[i] == NULL) {
                throw 102;
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
            sqltext = sqltext + std::to_string(parameters.experimental_profile[row][i]) + " ";
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
write_model_profile_to_db(sqlite3* db, double* D, double* ND, double* NS, double* species_4) //reading data using callback functions
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
                        INTENSITY_ARRAY = INTENSITY_ARRAY + std::to_string(D[x + row * parameters.X] / parameters.dye_conc[row]) + " ";
                    }
                    break;
                case 1:
                    species = "Bound Dye";
                    width = parameters.X;
                    for (int x = 0; x < width; x++) {
                        INTENSITY_ARRAY = INTENSITY_ARRAY + std::to_string(ND[x + row * parameters.X] / parameters.dye_conc[row]) + " ";
                    }
                    break;
                case 2:
                    species = "Double Bound Dye";
                    width = parameters.X;
                    for (int x = 0; x < width; x++) {
                        INTENSITY_ARRAY = INTENSITY_ARRAY + std::to_string(species_4[x + row * parameters.X] / parameters.dye_conc[row]) + " ";
                    }
                    break;
                case 3:
                    species = "Total Dye";
                    width = parameters.X;
                    for (int x = 0; x < width; x++) {
                        INTENSITY_ARRAY = INTENSITY_ARRAY + std::to_string((D[x + row * parameters.X] + ND[x + row * parameters.X] + species_4[x + row * parameters.X]) / parameters.dye_conc[row]) + " ";
                    }
                    break;
                case 4:
                    species = "Unbound Beads (wt%)";
                    width = parameters.X;
                    for (int x = 0; x < parameters.X; x++) {
                        INTENSITY_ARRAY = INTENSITY_ARRAY + std::to_string(NS[x + row * parameters.X] * 100.0d * 1.05d * ( 4.0d / 3.0d * M_PI * pow(parameters.diameter / 2.0d, 3.0d))) + " ";
                    }
                    break;
                case 5:
                    species = "Bound Beads (wt%)";
                    width = parameters.X;
                    for (int x = 0; x < width; x++) {
                        INTENSITY_ARRAY = INTENSITY_ARRAY + std::to_string(ND[x + row * parameters.X] / parameters.p1 * 100.0d * 1.05d * ( 4.0d / 3.0d * M_PI * pow(parameters.diameter / 2.0d, 3.0d))) + " ";
                    }
                    break;
                case 6:
                    species = "Double Bound Beads (wt%)";
                    width = parameters.X;
                    for (int x = 0; x < width; x++) {
                        INTENSITY_ARRAY = INTENSITY_ARRAY + std::to_string(species_4[x + row * parameters.X] / parameters.p1 / parameters.p2 * 100.0d * 1.05d * ( 4.0d / 3.0d * M_PI * pow(parameters.diameter / 2.0d, 3.0d))) + " ";
                    }
                    break; 
                case 7:
                    species = "Total Beads (wt%)";
                    width = parameters.X;
                    for (int x = 0; x < width; x++) {
                        INTENSITY_ARRAY = INTENSITY_ARRAY + std::to_string((NS[x + row * parameters.X] + ND[x + row * parameters.X] / parameters.p1 + species_4[x + row * parameters.X] / parameters.p1 / parameters.p2) * 100.0d * 1.05d * ( 4.0d / 3.0d * M_PI * pow(parameters.diameter / 2.0d, 3.0d)) ) + " ";
                    }
                    break;
                case 8:
                    species = "Experimental Derivative";
                    width = parameters.window_size;
                    for (int x = 0; x < (width - 1); x++) {
                        if ((x == 0) || (x + 1) == (width - 1)) {
                            INTENSITY_ARRAY = INTENSITY_ARRAY + "0 ";
                        } else {
                            INTENSITY_ARRAY = INTENSITY_ARRAY + std::to_string((parameters.experimental_profile[row][x + 1] - parameters.experimental_profile[row][x-1]) / (2 * parameters.W * 1.0e+6 / parameters.window_size)) + " ";
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
                            INTENSITY_ARRAY = INTENSITY_ARRAY + std::to_string((parameters.numeric_model_profile[row][x+1] / parameters.dye_conc[row] - parameters.numeric_model_profile[row][x - 1] / parameters.dye_conc[row]) / (2 * parameters.W * 1.0e+6 / parameters.window_size)) + " ";
                        }
                    }
                    break;    
                case 10:
                    species = "Experimental Profile";
                    width = parameters.window_size;
                    for (int x = 0; x < width; x++) {
                        INTENSITY_ARRAY = INTENSITY_ARRAY + std::to_string(parameters.experimental_profile[row][x]) + " ";
                    }
                    break;
                case 11:
                    species = "Numeric Model Profile";
                    width = parameters.window_size;
                    for (int x = 0; x < width; x++) {
                        INTENSITY_ARRAY = INTENSITY_ARRAY + std::to_string(parameters.numeric_model_profile[row][x] / parameters.dye_conc[row]) + " ";
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
normalize_profile() 
{
    double low_ref, high_ref, peak_ref;
    for (int j = 0; j < parameters.row_count; j++) {
        low_ref = 0.0d;
        high_ref = 0.0d;
        peak_ref = 0.0d;
        for (int i = parameters.low_ref_start + parameters.exp_left_padding; i < parameters.high_ref_end + parameters.exp_left_padding; i++) {
            if (i <= parameters.low_ref_end + parameters.exp_left_padding) {
                low_ref += parameters.experimental_profile[j][i];
            } else if (i >= parameters.high_ref_start + parameters.exp_left_padding) {
                high_ref += parameters.experimental_profile[j][i];
            }
        }
        for (int i = parameters.exp_left_padding; i < parameters.high_ref_end + parameters.exp_left_padding; i++) {
            peak_ref = max(parameters.experimental_profile[j][i], peak_ref);
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
                    parameters.experimental_profile[j][i] = 0;
                } else if (i >= parameters.high_ref_end + parameters.exp_left_padding) {  
                    parameters.experimental_profile[j][i] = 1;
                } else {
                    parameters.experimental_profile[j][i] = (parameters.experimental_profile[j][i] - low_ref) / (high_ref - low_ref);
                }
            } else if (parameters.normalization_method == "peak linear scaling") {
                if (i < parameters.low_ref_start + parameters.exp_left_padding) {
                    parameters.experimental_profile[j][i] = 0;
                } else if(i >= parameters.window_size - parameters.exp_right_padding) {
                    parameters.experimental_profile[j][i] = (high_ref - low_ref) / (peak_ref - low_ref);
                } else {
                    parameters.experimental_profile[j][i] = (parameters.experimental_profile[j][i] - low_ref) / (peak_ref - low_ref);
                }
            } else if (parameters.normalization_method == "none") {

            } 
        }
    }
}

double
scattering_correction(double NS, double ND, double p)
{
    if (parameters.scatter_correction_type == "NS_ND") {
        double bead_wt = (NS + ND / p) * 100.0d * 1.05d * ( 4.0d / 3.0d * M_PI * pow(parameters.diameter / 2.0d, 3.0d));
        if (parameters.diameter < 30.0e-9d) {
            return -0.238826108843563 * std::exp(-pow(0.0258645848310996 - bead_wt, 2.0d) / (2.0d * pow(0.00418992180425235, 2.0d))) + bead_wt * 1.40044738887211 + 1;
        };
        return -0.448191328804794 * std::exp(-pow(0.0711222856018783 - bead_wt, 2.0d) / (2.0d * pow(0.012623365952763, 2.0d))) + bead_wt * 2.14776259822044 + 1;
    } else {
        return 1.00d;
    }
}

void
solver(double* control_parameters, alglib::real_1d_array &residuals) // 
{
    
    //The 3 Concentration Arrays.
    double* D_in = parameters.D_in; 
    double* ND_in = parameters.ND_in; 
    double* NS_in = parameters.NS_in; 
    double* species_4_in = parameters.species_4_in;
    double* D = parameters.D; 
    double* NS = parameters.NS; 
    double* ND = parameters.ND; 
    double* species_4 = parameters.species_4;

    //initialize the 3 Concentration Arrays. 
    for (int row = 0; row < parameters.row_count; row++) {
        for (int z = 0; z < parameters.Z; z++) {
            for (int x = 0; x < parameters.X; x++) {
                int i = x + parameters.X * z;
                D[i] = D_in[x + row * parameters.X];
                NS[i] = NS_in[x + row * parameters.X];
                ND[i] = ND_in[x + row * parameters.X];
                species_4[i] = species_4_in[x + row * parameters.X];
            }
        }
    }
    double* D_out = parameters.D_out;
    double* NS_out = parameters.NS_out;
    double* ND_out = parameters.ND_out;
    double* species_4_out = parameters.species_4_out;
    double model_profile[parameters.row_count * parameters.window_size];

    double left_pad = parameters.left_pad, right_pad = parameters.right_pad;
    double p1 = parameters.p1, kon1 = parameters.kon1 * parameters.dt, koff1 = parameters.koff1 * parameters.dt, keq1 = parameters.keq1, QE1 = parameters.QE1;
    double p2 = parameters.p2, kon2 = parameters.kon2 * parameters.dt, koff2 = parameters.koff2 * parameters.dt, keq2 = parameters.keq2, QE2 = parameters.QE2;

    for (int index = 0; index < parameters.number_of_variables; index++) { // This sorts out out of order parameters
        if (parameters.solve_for[index] == "left_pad") {
            left_pad    = control_parameters[index];
        } else if (parameters.solve_for[index] == "right_pad") {
            right_pad   = control_parameters[index];
        } else if (parameters.solve_for[index] == "kon1") {
            kon1        = control_parameters[index] * parameters.dt;
        } else if (parameters.solve_for[index] == "p1") {
            p1          = control_parameters[index];
        } else if (parameters.solve_for[index] == "koff1") {
            koff1       = control_parameters[index] * parameters.dt;
        } else if (parameters.solve_for[index] == "keq1") {
            keq1        = control_parameters[index];
        } else if (parameters.solve_for[index] == "QE1") {
            QE1         = control_parameters[index];
        } else if (parameters.solve_for[index] == "kon2") {
            kon2        = control_parameters[index] * parameters.dt;
        } else if (parameters.solve_for[index] == "p2") {
            p2          = control_parameters[index];
        } else if (parameters.solve_for[index] == "koff2") {
            koff2       = control_parameters[index] * parameters.dt;
        } else if (parameters.solve_for[index] == "keq2") {
            keq2        = control_parameters[index];
        } else if (parameters.solve_for[index] == "QE2") {
            QE2         = control_parameters[index];
        }
    }
    if (!parameters.solve_for_koff1) {koff1 = kon1 / keq1;}
    if (!parameters.solve_for_koff2) {koff2 = kon2 / keq2;}
    
    double scatter = 1.0d;

    Eigen::MatrixXd m(parameters.X,parameters.X);
    Eigen::MatrixXd m_beads(parameters.X,parameters.X);
    Eigen::MatrixXd E_D(parameters.X,1);
    Eigen::MatrixXd E_NS(parameters.X,1);
    Eigen::MatrixXd E_ND(parameters.X,1);
    Eigen::MatrixXd E_species_4(parameters.X,1);
    Eigen::MatrixXd solution_D(parameters.X,1);
    Eigen::MatrixXd solution_NS(parameters.X,1);
    Eigen::MatrixXd solution_ND(parameters.X,1);
    Eigen::MatrixXd solution_species_4(parameters.X,1);

    for(int j = 0; j < parameters.X; j++) {
        for (int i = 0; i < parameters.X; i++) {
            m(j,i)              = (*matrixes.m)(j,i);
            m_beads(j,i)        = (*matrixes.m_beads)(j,i);
        }
        E_D(j)                  = (*matrixes.E_D)(j);
        E_NS(j)                 = (*matrixes.E_NS)(j);
        E_ND(j)                 = (*matrixes.E_ND)(j);
        E_species_4(j)          = (*matrixes.E_species_4)(j);
        solution_D(j)           = (*matrixes.solution_D)(j);
        solution_NS(j)          = (*matrixes.solution_NS)(j);
        solution_ND(j)          = (*matrixes.solution_ND)(j);
        solution_species_4(j)   = (*matrixes.solution_species_4)(j);
    }

    double reaction_rate1 = 0.0d, reaction_rate2 = 0.0d;
    for (int row = 0; row < parameters.row_count; row++) {
        int current_point = 0;
        for (int z = 0; z < parameters.Z; z++) {
            for (int i = 0; i < parameters.X; i++) {
                current_point = i + z * parameters.X;
                
                reaction_rate1 = kon1 * (D[current_point] * NS[current_point]) - koff1 * ND[current_point];
                reaction_rate2 = kon2 * (D[current_point] * ND[current_point]) - koff2 * species_4[current_point];
                reaction_rate1 = min(reaction_rate1, NS[current_point] * p1);
                reaction_rate2 = max(-species_4[current_point], reaction_rate2);
            
                if (-reaction_rate1 + reaction_rate2 / p2 > ND[current_point]) {
                    if (reaction_rate1 < 0 && reaction_rate2 > 0) {
                        reaction_rate1 = ND[current_point] * (-reaction_rate1) / (-reaction_rate1 + reaction_rate2 / p2);
                        reaction_rate2 = ND[current_point] * p2 * (reaction_rate2 / p2) / (-reaction_rate1 + reaction_rate2 / p2);
                    } else if (reaction_rate1 < 0) {
                        reaction_rate1 = -ND[current_point] - reaction_rate2 / p2;
                    } else if (reaction_rate2 > 0) {
                        reaction_rate2 = ND[current_point] * p2 + reaction_rate1;
                    }
                }
                if (reaction_rate1 + reaction_rate2 > D[current_point]) {
                    if (reaction_rate1 > 0 && reaction_rate2 > 0) {
                        reaction_rate1 = D[current_point] * reaction_rate1 / (reaction_rate1 + reaction_rate2);
                        reaction_rate2 = D[current_point] * reaction_rate2 / (reaction_rate1 + reaction_rate2);
                    } else if (reaction_rate1 > 0) {
                        reaction_rate1 = D[current_point] - reaction_rate2;
                    } else if (reaction_rate2 > 0) {
                        reaction_rate2 = D[current_point] - reaction_rate1;
                    }
                }
                
                if (parameters.disable_reactions) {
                    reaction_rate1 = 0.0d;
                    reaction_rate2 = 0.0d;
                }

                if (i == 0) { // left edge
                    E_D(i,0)            =                                                        (1.00d - parameters.r_dye)   * D[current_point]                  + parameters.r_dye   * D[current_point + 1]            - reaction_rate1 - reaction_rate2;
                    E_NS(i,0)           =                                                        (1.00d - parameters.r_beads) * NS[current_point]                 + parameters.r_beads * NS[current_point + 1]           - reaction_rate1 / p1;
                    E_ND(i,0)           =                                                        (1.00d - parameters.r_beads) * ND[current_point]                 + parameters.r_beads * ND[current_point + 1]           + reaction_rate1 - reaction_rate2 / p2;
                    E_species_4(i,0)    =                                                        (1.00d - parameters.r_beads) * species_4[current_point]          + parameters.r_beads * species_4[current_point + 1]    + reaction_rate2;
                } else if (i == parameters.X - 1) { // right edge
                    E_D(i,0)            = parameters.r_dye   * D[current_point - 1]            + (1.00d - parameters.r_dye)   * D[current_point]                                                                         - reaction_rate1 - reaction_rate2;
                    E_NS(i,0)           = parameters.r_beads * NS[current_point - 1]           + (1.00d - parameters.r_beads) * NS[current_point]                                                                        - reaction_rate1 / p1;
                    E_ND(i,0)           = parameters.r_beads * ND[current_point - 1]           + (1.00d - parameters.r_beads) * ND[current_point]                                                                        + reaction_rate1 - reaction_rate2 / p2;
                    E_species_4(i,0)    = parameters.r_beads * species_4[current_point - 1]    + (1.00d - parameters.r_beads) * species_4[current_point]                                                                 + reaction_rate2;
                } else { // mid points
                    E_D(i,0)            = parameters.r_dye   * D[current_point - 1]            + (2.00d - 2.00d * parameters.r_dye)   * D[current_point]          + parameters.r_dye   * D[current_point + 1]            - reaction_rate1 - reaction_rate2;
                    E_NS(i,0)           = parameters.r_beads * NS[current_point - 1]           + (2.00d - 2.00d * parameters.r_beads) * NS[current_point]         + parameters.r_beads * NS[current_point + 1]           - reaction_rate1 / p1;
                    E_ND(i,0)           = parameters.r_beads * ND[current_point - 1]           + (2.00d - 2.00d * parameters.r_beads) * ND[current_point]         + parameters.r_beads * ND[current_point + 1]           + reaction_rate1 - reaction_rate2 / p2;
                    E_species_4(i,0)    = parameters.r_beads * species_4[current_point - 1]    + (2.00d - 2.00d * parameters.r_beads) * species_4[current_point]  + parameters.r_beads * species_4[current_point + 1]    + reaction_rate2;
                }
            }

            solution_D              = m * E_D;
            solution_NS             = m_beads * E_NS;
            solution_ND             = m_beads * E_ND;
            solution_species_4      = m_beads * E_species_4;

            for (int i = 0; i < parameters.X; i++) {
                current_point = i + z * parameters.X;

                D[current_point]                          = max(min(solution_D(i,0), parameters.dye_conc[row]), 0.00d);
                D_out[i + row * parameters.X]             = parameters.dye_conc[row]; //D[current_point];
                NS[current_point]                         = max(min(solution_NS(i,0), parameters.bead_conc[row]), 0.00d);
                NS_out[i + row * parameters.X]            = parameters.bead_conc[row]; //NS[current_point];
                ND[current_point]                         = max(solution_ND(i,0), 0.00d);
                ND_out[i + row * parameters.X]            = ND[current_point];
                species_4[current_point]                  = max(solution_species_4(i,0), 0.00d);
                species_4_out[i + row * parameters.X]     = species_4[current_point];
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
                    scatter = scattering_correction(NS[bottom_point + (parameters.Z - 1) * parameters.X], ND[bottom_point + (parameters.Z - 1) * parameters.X], p1);
                    model_profile[row * parameters.window_size + i] = (D[bottom_point + (parameters.Z - 1) * parameters.X] + ND[bottom_point + (parameters.Z - 1) * parameters.X] * QE1 + species_4[bottom_point + (parameters.Z - 1) * parameters.X] * QE2) * scatter;
                } else {
                    scatter = scattering_correction(NS[bottom_point + (parameters.Z - 1) * parameters.X] * ((double)top_point - split) + NS[top_point + (parameters.Z - 1) * parameters.X] * (split - (double)bottom_point), ND[bottom_point + (parameters.Z - 1) * parameters.X] * ((double)top_point - split) + ND[top_point + (parameters.Z - 1) * parameters.X] * (split - (double)bottom_point), p1);
                    model_profile[row * parameters.window_size + i] = ((D[bottom_point + (parameters.Z - 1) * parameters.X] + ND[bottom_point + (parameters.Z - 1) * parameters.X] * QE1 + species_4[bottom_point + (parameters.Z - 1) * parameters.X] * QE2) * ((double)top_point - split) + (D[top_point + (parameters.Z - 1) * parameters.X] + ND[top_point + (parameters.Z - 1) * parameters.X] * QE1 + species_4[top_point + (parameters.Z - 1) * parameters.X] * QE2) * (split - (double)bottom_point)) * scatter;
                }
            } else {
                model_profile[row * parameters.window_size + i] = parameters.dye_conc[row];
            }
            residuals[i + row * parameters.window_size] = pow(model_profile[row * parameters.window_size + i] / parameters.dye_conc[row] - parameters.experimental_profile[row][i], 2.0d);
            parameters.numeric_model_profile[row][i] = model_profile[row * parameters.window_size + i];
        }
    }
}

void
alglib_solver(const alglib::real_1d_array &control_parameters, alglib::real_1d_array &residuals, void *ptr)
{
    parameters.iterations = parameters.iterations + 1;
    double control_parameters_temp [parameters.number_of_variables];
    for (int i = 0; i < parameters.number_of_variables; i++) {
        control_parameters_temp [i] = control_parameters[i];
    }
    solver(control_parameters_temp, residuals);
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
    std::cout << " Done" << endl;
    parameters.file_name = "../" + parameters.experiment_name;
    parameters.M = parameters.X * parameters.Z;

    std::cout << "Retrieving line count:";
    lines_from_profile_text(db);
    std::cout << " Done" << endl;

    parameters.wt_percent                = new double[parameters.row_count]();   //= 0.1d; // wt%
    parameters.dye_conc_mgml             = new double[parameters.row_count]();   //= 0.00336d; // mg/ml FITC
    parameters.bead_conc                 = new double[parameters.row_count]();   //= wt_percent / 100.0d / 1.05d / ( 4.0d / 3.0d * M_PI * pow(diameter / 2.0d, 3.0d)); // beads/m3 
    parameters.dye_conc                  = new double[parameters.row_count]();   //= dye_conc_mgml * 1000.0d / 332.326d * 6.022e+23;     // molecules FITC / m3
    parameters.beginning_of_channel      = new int[parameters.row_count]();
    parameters.end_of_channel            = new int[parameters.row_count]();

    parameters.ENTRANCE_1_SPECIES_1_CONC = new double[parameters.row_count](); 
    parameters.ENTRANCE_2_SPECIES_1_CONC = new double[parameters.row_count](); 
    parameters.ENTRANCE_1_SPECIES_2_CONC = new double[parameters.row_count](); 
    parameters.ENTRANCE_2_SPECIES_2_CONC = new double[parameters.row_count]();

    std::cout << "Reading experimental parameters from database: ";
    read_exp_parameters_from_db(db);
    std::cout << "Done" << endl;
    parameters.output_file_name = parameters.file_name + " " + parameters.normalization_method + " scatter_" + parameters.scatter_correction_type + " output.txt";

    parameters.scale                     = new double[parameters.number_of_variables];
    parameters.initial_values_alglib     = new double[parameters.number_of_variables];
    parameters.low_bound                 = new double[parameters.number_of_variables];
    parameters.up_bound                  = new double[parameters.number_of_variables];

    std::cout << "Reading experimental profiles from database: ";
    read_raw_profiles_from_db(db);
    std::cout << "Done" << endl;

    std::cout << "Reading alglib inputs from database: ";
    read_alglib_values_from_db(db);
    std::cout << "Done" << endl;

    std::cout << "scale:";
    for (int item = 0; item < parameters.number_of_variables; item++) {
        std::cout << parameters.scale[item] << " ";
    }
    std::cout << endl;

    std::cout << "Normalizing " << parameters.row_count << " profiles with " << parameters.normalization_method << ": ";
    normalize_profile();
    std::cout << "Done" << endl;

    std::cout << "Writing " << parameters.row_count << " normalized profiles to database: ";
    write_normalized_values_to_db(db);
    std::cout << "Done" << endl;

    std::cout << "Priming " << parameters.output_file_name << ": ";
    prime_excel_output(parameters.output_file_name);
    std::cout << "Done" << endl;

    parameters.visc            = 0.0010016d; // Dynamic viscosity of water at 20C in Pa.s
    parameters.difusion_dye    = 4.9e-10d; // m2/s From 4.9 × 10−6 cm2 s−1 The diffusion coefficient of fluorescein in water at 21.5°C, as calculated from the Wilke-Chang correlation
    parameters.difusion_beads  = 1.380649e-23d * (273.15d + 25.0d) / (3.0d * M_PI * parameters.visc * parameters.diameter); // m2/s From kB*T/(3*pi*visc*d) kB=1.380649×10−23 J⋅K−1

    // Derived values for Crank-Nicolson implicit method
    parameters.restime       = parameters.W * parameters.H * parameters.L / parameters.flow; // seconds
    parameters.dt            = parameters.restime / parameters.Z; // seconds
    parameters.T_dye         = parameters.difusion_dye * parameters.restime / (parameters.W * parameters.W);
    parameters.T_beads       = parameters.difusion_beads * parameters.restime / (parameters.W * parameters.W);
    parameters.dT_dye        = parameters.difusion_dye * parameters.dt / (parameters.W * parameters.W); // T = Dt/l^2
    parameters.dT_beads      = parameters.difusion_beads * parameters.dt / (parameters.W * parameters.W);
    parameters.dx            = parameters.W / parameters.X; // m / subdivision
    parameters.dX            = 1.0d / parameters.X; // X = x/l
    parameters.r_dye         = parameters.dT_dye / parameters.dX / parameters.dX; // r = dT/(dX)^2
    parameters.r_beads       = parameters.dT_beads / parameters.dX / parameters.dX; // r = dT/(dX)^2

    //The 3 Concentration Arrays.
    double* D_in{new double[parameters.row_count * parameters.X]{}};
    double* ND_in{new double[parameters.row_count * parameters.X]{}};
    double* NS_in{new double[parameters.row_count * parameters.X]{}};
    double* species_4_in{new double[parameters.row_count * parameters.X]{}};
    parameters.D_in = D_in;
    parameters.NS_in = NS_in;
    parameters.ND_in = ND_in;
    parameters.species_4_in = species_4_in;
    double* D{new double[parameters.M]{}};
    double* ND{new double[parameters.M]{}};
    double* NS{new double[parameters.M]{}};
    double* species_4{new double[parameters.M]{}};
    parameters.D = D;
    parameters.ND = ND;
    parameters.NS = NS;
    parameters.species_4 = species_4;
    double* D_out{new double[parameters.X * parameters.row_count]{}};
    double* ND_out{new double[parameters.X * parameters.row_count]{}};
    double* NS_out{new double[parameters.X * parameters.row_count]{}};
    double* species_4_out{new double[parameters.X * parameters.row_count]{}};
    parameters.D_out = D_out;
    parameters.ND_out = ND_out;
    parameters.NS_out = NS_out;
    parameters.species_4_out = species_4_out;

    for (int row = 0; row < parameters.row_count; row++) {
        //initialize the concentration inlet arrays. 
        for (int x = 0; x < parameters.X; x++) {
            int i = x + parameters.X * row;
            if (x < (static_cast<int>(parameters.X * parameters.ENTRANCE_1_FLOWRATE / (parameters.ENTRANCE_1_FLOWRATE + parameters.ENTRANCE_2_FLOWRATE)))) {
                if (parameters.ENTRANCE_1_SPECIES_1_NAME == "FITC") {
                    D_in[i] = parameters.ENTRANCE_1_SPECIES_1_CONC[row] * 1000.0d / 332.326d * 6.022e+23;  // molecules / m3 
                } else if (parameters.ENTRANCE_1_SPECIES_2_NAME == "FITC") {
                    D_in[i] = parameters.ENTRANCE_1_SPECIES_2_CONC[row] * 1000.0d / 332.326d * 6.022e+23;  // molecules / m3
                } else {
                    D_in[i]  = 0.0d;
                }
                if (parameters.ENTRANCE_1_SPECIES_1_NAME == "20nm PS" || parameters.ENTRANCE_1_SPECIES_1_NAME == "40nm PS" || parameters.ENTRANCE_1_SPECIES_2_NAME == "20nm PS" || parameters.ENTRANCE_1_SPECIES_2_NAME == "40nm PS") {
                    NS_in[i] = parameters.bead_conc[row];
                } else {
                    NS_in[i] = 0.0d;
                }
            } else {
                if (parameters.ENTRANCE_2_SPECIES_1_NAME == "FITC") {
                    D_in[i] = parameters.ENTRANCE_2_SPECIES_1_CONC[row] * 1000.0d / 332.326d * 6.022e+23;  // molecules / m3 
                } else if (parameters.ENTRANCE_2_SPECIES_2_NAME == "FITC") {
                    D_in[i] = parameters.ENTRANCE_2_SPECIES_2_CONC[row] * 1000.0d / 332.326d * 6.022e+23;  // molecules / m3 
                } else {
                    D_in[i]  = 0.0d;
                }
                if (parameters.ENTRANCE_2_SPECIES_1_NAME == "20nm PS" || parameters.ENTRANCE_2_SPECIES_1_NAME == "40nm PS" || parameters.ENTRANCE_2_SPECIES_2_NAME == "20nm PS" || parameters.ENTRANCE_2_SPECIES_2_NAME == "40nm PS") {
                    NS_in[i] = parameters.bead_conc[row];
                } else {
                    NS_in[i] = 0.0d;
                }
            }
            ND_in[i] = 0.0d;
            species_4_in[i] = 0.0d;

        }
    }
    Eigen::MatrixXd m                   = Eigen::MatrixXd::Zero(parameters.X, parameters.X);
    matrixes.m                          = &m;
    Eigen::MatrixXd m_beads             = Eigen::MatrixXd::Zero(parameters.X, parameters.X);
    matrixes.m_beads                    = &m_beads;
    Eigen::MatrixXd E_D                 = Eigen::MatrixXd::Zero(1, parameters.X);
    matrixes.E_D                        = &E_D;
    Eigen::MatrixXd E_NS                = Eigen::MatrixXd::Zero(1, parameters.X);
    matrixes.E_NS                       = &E_NS;
    Eigen::MatrixXd E_ND                = Eigen::MatrixXd::Zero(1, parameters.X);
    matrixes.E_ND                       = &E_ND;
    Eigen::MatrixXd E_species_4         = Eigen::MatrixXd::Zero(1, parameters.X);
    matrixes.E_species_4                = &E_species_4;
    Eigen::MatrixXd solution_D          = Eigen::MatrixXd::Zero(1, parameters.X);
    matrixes.solution_D                 = &solution_D;
    Eigen::MatrixXd solution_NS         = Eigen::MatrixXd::Zero(1, parameters.X);
    matrixes.solution_NS                = &solution_NS;
    Eigen::MatrixXd solution_ND         = Eigen::MatrixXd::Zero(1, parameters.X);
    matrixes.solution_ND                = &solution_ND;
    Eigen::MatrixXd solution_species_4  = Eigen::MatrixXd::Zero(1, parameters.X);
    matrixes.solution_species_4         = &solution_species_4;

    for (int i = 0; i < parameters.X; i++) {
        for (int j = 0; j < parameters.X; j++) {
            if (j == i - 1) {
                m(i , j) = -parameters.r_dye;
                m_beads(i , j) = -parameters.r_beads;
            } else if (j == i + 1) {
                m(i , j) = -parameters.r_dye;
                m_beads(i , j) = -parameters.r_beads;
            } else if ((i != 0 && i != parameters.X - 1) && (i == j)) {
                m(i , j) = 2.0d + 2.0d * parameters.r_dye;
                m_beads(i , j) = 2.0d + 2.0d * parameters.r_beads;
            } else if ((i == 0 || i == parameters.X - 1) && (i == j)) {
                m(i , j) = 1.0d + parameters.r_dye;
                m_beads(i , j) = 1.0d + parameters.r_beads;
            } else {
                m(i , j) = 0.0d;
                m_beads(i , j) = 0.0d;
            }
        }
    }
    m = m.inverse();
    m_beads = m_beads.inverse();

    for (int i = 0; i < parameters.X * parameters.row_count; i++) {
        D_out[i] = 0;
        ND_out[i] = 0;
        NS_out[i] = 0;
        species_4_out[i] = 0;
    }

    try
    {
        double epsx = 1e-9;
        /*
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
        alglib::ae_int_t maxits = 0;
        
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
        alglib::minlmsetcond(state, epsx, maxits);
        alglib::minlmsetscale(state, s);
        alglib::minlmsetnonmonotonicsteps(state, 2);
        if (parameters.run_solver) {
            std::cout << "minlmoptimize: " << parameters.experiment_name << " ";
            alglib::minlmoptimize(state, alglib_solver);   // Optimize
            std::cout << "Done after " << parameters.iterations << " iterations" << endl;
            alglib::minlmresults(state, control_parameters, rep);
            printf("%s\n", control_parameters.tostring(4).c_str());
        };

        for (int index = 0; index < parameters.number_of_variables; index++) { // This sorts out of order parameters
            if (parameters.solve_for[index] == "left_pad") {
                parameters.left_pad    = control_parameters[index];
            } else if (parameters.solve_for[index] == "right_pad") {
                parameters.right_pad   = control_parameters[index];
            } else if (parameters.solve_for[index] == "p1") {
                parameters.p1          = control_parameters[index];
            } else if (parameters.solve_for[index] == "kon1") {
                parameters.kon1        = control_parameters[index];
            } else if (parameters.solve_for[index] == "koff1") {
                parameters.koff1       = control_parameters[index];
            } else if (parameters.solve_for[index] == "keq1") {
                parameters.keq1        = control_parameters[index];
            } else if (parameters.solve_for[index] == "QE1") {
                parameters.QE1         = control_parameters[index];
            } else if (parameters.solve_for[index] == "p2") {
                parameters.p2          = control_parameters[index];
            } else if (parameters.solve_for[index] == "kon2") {
                parameters.kon2        = control_parameters[index];
            } else if (parameters.solve_for[index] == "koff2") {
                parameters.koff2       = control_parameters[index];
            } else if (parameters.solve_for[index] == "keq2") {
                parameters.keq2        = control_parameters[index];
            } else if (parameters.solve_for[index] == "QE2") {
                parameters.QE2         = control_parameters[index];
            }
        }
        if (parameters.solve_for_keq1) {parameters.koff1 = parameters.kon1 / parameters.keq1;}
        if (parameters.solve_for_keq2) {parameters.koff2 = parameters.kon2 / parameters.keq2;}
 
        save_excel_output(D_out, ND_out, NS_out, species_4_out);
        std::cout << "Writing model profiles to database: ";
        write_model_profile_to_db(db, D_out, ND_out, NS_out, species_4_out);
        std::cout << "Done" << endl;
    }
    catch(alglib::ap_error alglib_exception)
    {
        printf("ALGLIB exception with message '%s'\n", alglib_exception.msg.c_str());
        return 1;
    }

    sqlite3_close(db);
    fprintf(stdout, "Database closed.\n");

    std::cout << "Program finished " << endl;
    return 0;
}