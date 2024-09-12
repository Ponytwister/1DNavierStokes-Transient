#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include <eigen-3.4.0/Eigen/Dense>
#include <optimization.h>
#include <stdafx.h>
#include <solvers.h>
// D: Dye
// NS: Nanoparticle Site
// ND: Nanoparticle Bound Dye

using namespace std;

void prime_excel_output(string file_name);
void save_excel_output(double* D, double* ND, double* NS, double* control_parameters);
int lines_from_profile_text(string file_name);
int max_profile_array_size_int(int line, string file_name);
void read_profile_from_text(string file_name);
void normalize_profile();
void solver(double* control_parameters, alglib::real_1d_array &residuals); 
void alglib_shift_solver(const alglib::real_1d_array &control_parameters, alglib::real_1d_array &residuals, void *ptr);
void alglib_react_solver(const alglib::real_1d_array &control_parameters, alglib::real_1d_array &residuals, void *ptr);
void alglib_solver(const alglib::real_1d_array &control_parameters, alglib::real_1d_array &residuals, void *ptr);

const int Z                 = 100; // divisions along the z(time) axis
const int X                 = 500; // divisions along the x axis
const int M                 = X * Z;

// Device Dimenssions
const double W = 5e-4;  //meters: 500 um
const double H = 4e-5;  //meters: 40 um
const double L = 0.025; //meters: 2.5 cm

// Operating conditions/settings
const double flow            = 2.0d * 5e-9d / 60;                        //m3/s: 2*5 ulmin
const double diameter        = 23e-9d; // meters
string file_name = "../../ExampleProfiles20nm24_3rdNormalized";
string input_file_name = file_name + ".txt";
string output_file_name = file_name + "_output.txt";
const int profile_file_line_count = lines_from_profile_text(input_file_name);
const int window_size = 180;
const int exp_left_padding = 20;
const int exp_right_padding = 20;

double num_left_padding = 39.7704d;
double num_right_padding = 20.0155d;
double num_profile_width = window_size - num_left_padding - num_right_padding;

int low_ref_start = 0;
int low_ref_end = 20;
int high_ref_start = 130;
int high_ref_end = 140;

// Physical Constants
const double visc            = 0.0010016d; // Dynamic viscosity of water at 20C in Pa.s
const double difusion_dye    = 4.9e-10d; // m2/s From 4.9 × 10−6 cm2 s−1 The diffusion coefficient of fluorescein in water at 21.5°C, as calculated from the Wilke-Chang correlation
const double difusion_beads  = 1.380649e-23d * (273.15d + 25.0d) / (3.0d * M_PI * visc * diameter); // m2/s From kB*T/(3*pi*visc*d) kB=1.380649×10−23 J⋅K−1

// Derived values for Crank-Nicolson implicit method
const double restime       = W * H * L / flow; // seconds
const double dt            = restime / Z; // seconds
const double T_dye         = difusion_dye * restime / (W * W);
const double T_beads       = difusion_beads * restime / (W * W);
const double dT_dye        = difusion_dye * dt / (W * W); // T = Dt/l^2
const double dT_beads      = difusion_beads * dt / (W * W);
const double dx            = W / X; // m / subdivision
const double dX            = 1.0d / X; // X = x/l
const double r_dye         = dT_dye / dX / dX; // r = dT/(dX)^2
const double r_beads       = dT_beads / dX / dX; // r = dT/(dX)^2

struct parameters_struct {
    double p               = 124.2129d; // FITC molecules / PS bead
    double kon             = 1.987E-14d * (W / X); // 0.10d;
    double koff            = 1.987E-16d * (W / X); 
    double* wt_percent;                      //= 0.1d; // wt%
    double* dye_conc_mgml;                   //= 0.00336d; // mg/ml FITC
    double* bead_conc;                       //= wt_percent / 100.0d / 1.05d / ( 4.0d / 3.0d * M_PI * powf(diameter / 2.0d, 3.0d)); // beads/m3 
    double* dye_conc;                        //= dye_conc_mgml * 1000.0d / 332.326d * 6.022e+23;     // molecules FITC / m3
    int* beginning_of_channel;              //= 56;
    int* end_of_channel;                    //= 196;
    int* profile_array_size_int;            //= 140;
    int* max_profile_array_size_int;
    double experimental_profile[16][180];    // 15_111 16_140
    double numeric_model_profile[16][180];
} parameters;

struct matrix_struct {
    Eigen::Matrix<double, -1, -1>* m;
    Eigen::Matrix<double, -1, -1>* m_beads;
    Eigen::Matrix<double, -1, -1>* E;
    Eigen::Matrix<double, -1, -1>* E_NS;
    Eigen::Matrix<double, -1, -1>* E_ND;
    Eigen::Matrix<double, -1, -1>* solution;
    Eigen::Matrix<double, -1, -1>* solution_NS;
    Eigen::Matrix<double, -1, -1>* solution_ND;
} matrixes;

void
prime_excel_output(string file_name)
{
    ofstream fout;
    fout.open(file_name, std::ofstream::out | std::ofstream::trunc);
    fout << "res_time  " << "bind_ratio(p)  " << "forward_reaction_rate  "  << "reverse_reaction_rate  "<< "dye_conc.   " << "bead_conc.  "<< "species  "; 
    for (int x = 0; x < X; x++) {
        fout << x << "  ";
    }
    fout << endl;
    fout << "sec  " << "molc/bead  " << "forward_reaction_rate  " << "reverse_reaction_rate  " << "molc/m3   " << "bead/m3  "<< "Channel_Width_(um)->  "; 
}

void
save_excel_output(double* D, double* ND, double* NS, double* control_parameters)
{
    double dye_bead_ratio = 2000.0d;
    double total_dye = 0.0d;
    int z = Z - 1;
    ofstream fout;
    fout.open(output_file_name, std::ofstream::out | std::ofstream::app);

    double scale_factor = W * 1.0e6 / ((double)window_size - control_parameters[1] - control_parameters[0]);
    for (int x = 0; x < X; x++) {
        fout << (x  - control_parameters[0]) * scale_factor << "  ";
    }
    fout << endl;

    for (int step_in_line_read = 0; step_in_line_read < profile_file_line_count; step_in_line_read++) {
        for (int j = 0; j < 9; j++) {
            fout << z * dt << "  " << control_parameters[2] << "  " << control_parameters[3] << "  " << control_parameters[4] << "  " << parameters.dye_conc[step_in_line_read] << "   " <<  parameters.bead_conc[step_in_line_read] << "  ";
            switch(j) {
                case 0:
                    fout << "Free_Dye  ";
                    for (int x = 0; x < X; x++) {
                        fout << D[x + step_in_line_read * X] << "    ";
                    }
                    break;
                case 1:
                    fout << "Bound_Dye  ";
                    for (int x = 0; x < X; x++) {
                        fout << ND[x + step_in_line_read * X] << "    ";
                    }
                    break;
                case 2:
                    fout << "Total_Dye  ";
                    for (int x = 0; x < X; x++) {
                        fout << D[x + step_in_line_read * X] + ND[x + step_in_line_read * X] << "    ";
                    }
                    break;
                case 3:
                    fout << "Unbound_Beads  ";
                    for (int x = 0; x < X; x++) {
                        fout << NS[x + step_in_line_read * X] << "    ";
                    }
                    break;
                case 4:
                    fout << "Bound_Beads  ";
                    for (int x = 0; x < X; x++) {
                        fout << ND[x + step_in_line_read * X] / control_parameters[2] << "    ";
                    }
                    break;   
                case 5:
                    fout << "Total_Beads  ";
                    for (int x = 0; x < X; x++) {
                        fout << NS[x + step_in_line_read * X] + ND[x + step_in_line_read * X] / control_parameters[2] << "    ";
                    }
                    break;
                case 6:
                    fout << "Exitation  ";
                    for (int x = 0; x < X; x++) {
                        dye_bead_ratio = (D[x + step_in_line_read * X] + ND[x + step_in_line_read * X]) / (NS[x+ step_in_line_read * X] + ND[x + step_in_line_read * X] / parameters.p);
                        total_dye = D[x + step_in_line_read * X] + ND[x + step_in_line_read * X];
                        if(isinf(dye_bead_ratio)) {
                            fout << total_dye / parameters.dye_conc[step_in_line_read] << "    ";
                        }
                        else if(dye_bead_ratio < 1000) {
                            fout << total_dye * (-0.00137 * dye_bead_ratio + 2.04) / parameters.dye_conc[step_in_line_read] << "    ";
                        }
                        else if(dye_bead_ratio < 1700) {
                            fout << max(min(total_dye * (-0.000000359 * powf(dye_bead_ratio, 2.0d) + 0.00147 * dye_bead_ratio - 0.432d), 1.0) , total_dye * (-0.00137 * dye_bead_ratio + 2.04)) / parameters.dye_conc[step_in_line_read] << "    ";
                        }
                        else {
                            fout << total_dye / parameters.dye_conc[step_in_line_read] << "    ";
                        }
                    }
                    break;
                case 7:
                    fout << "Experimental_Profile  ";
                    for (int x = 0; x <  parameters.max_profile_array_size_int[step_in_line_read]; x++) {
                        fout << parameters.experimental_profile[step_in_line_read][x] << "    ";
                    }
                    break;
                case 8:
                    fout << "Total_Dye_rescale  ";
                    for (int x = 0; x <  parameters.max_profile_array_size_int[step_in_line_read]; x++) {
                        fout << parameters.numeric_model_profile[step_in_line_read][x] / parameters.dye_conc[step_in_line_read] << "    ";
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
lines_from_profile_text(string file_name)
{
    fstream fin;
    fin.open(file_name, std::fstream::in | std::fstream::binary); 
    string data;
    int line_count = 0;
    fin.seekg(0);
    while (!fin.eof()) { 
        fin >> data;
        if(data == "end") {
            line_count++;
        }
    }
    return line_count - 1;
}

int
max_profile_array_size_int(string file_name)
{
    fstream fin;
    fin.open(file_name, std::fstream::in | std::fstream::binary); 
    string data;
    int array_size = 0;
    fin.seekg(0); 
    fin >> data;
    fin >> data;
    fin >> data;
    array_size = -stoi(data);
    fin >> data;
    array_size = array_size + stoi(data);
    return array_size;
}

void
read_profile_from_text(string file_name)
{
    fstream fin;
    fin.open(file_name, std::fstream::in | std::fstream::binary); 
    string data;
    int step_in_line_read = 0;
    fin.seekg(0);
    fin >> data;
    low_ref_start = stoi(data);
    fin >> data;
    low_ref_end = stoi(data);
    fin >> data;
    high_ref_start = stoi(data);
    fin >> data;
    high_ref_end = stoi(data);
    fin >> data;
    while (!fin.eof()) {    
        parameters.wt_percent[step_in_line_read] = stod(data);
        parameters.bead_conc[step_in_line_read] = parameters.wt_percent[step_in_line_read] / 100.0d / 1.05d / ( 4.0d / 3.0d * M_PI * powf(diameter / 2.0d, 3.0d)); // beads/m3 

        fin >> data;
        parameters.dye_conc_mgml[step_in_line_read] = stod(data);
        parameters.dye_conc[step_in_line_read] = parameters.dye_conc_mgml[step_in_line_read] * 1000.0d / 332.326d * 6.022e+23;     // molecules FITC / m3
        
        fin >> data;
        parameters.beginning_of_channel[step_in_line_read] = stoi(data);
        
        fin >> data;
        parameters.end_of_channel[step_in_line_read] = stoi(data);
        parameters.profile_array_size_int[step_in_line_read] = parameters.end_of_channel[step_in_line_read] - parameters.beginning_of_channel[step_in_line_read];
        parameters.max_profile_array_size_int[step_in_line_read] = parameters.profile_array_size_int[step_in_line_read] + exp_left_padding + exp_right_padding;
        fin >> data;
        for (int i=0; data != "end"; i++) {
            if(i >= parameters.beginning_of_channel[step_in_line_read] && i < parameters.end_of_channel[step_in_line_read]) {
                parameters.experimental_profile[step_in_line_read][exp_left_padding + i - parameters.beginning_of_channel[step_in_line_read]] = stod(data);
            }
            fin >> data;
        }
        if(!fin.eof()) {
            fin >> data;
        }
        step_in_line_read++;
    }
}

void 
normalize_profile() 
{
    double low_ref;
    double high_ref;
    for (int j = 0; j < profile_file_line_count; j++) {
        low_ref = 0.0d;
        high_ref = 0.0d;
        for (int i = low_ref_start + exp_left_padding; i < high_ref_end + exp_left_padding; i++) {
            if (i <= low_ref_end + exp_left_padding) {
                low_ref += parameters.experimental_profile[j][i];
            } else if (i >= high_ref_start + exp_left_padding) {
                high_ref += parameters.experimental_profile[j][i];
            }
        }
        low_ref = low_ref / (low_ref_end - low_ref_start + 1);
        high_ref = high_ref / (high_ref_end - high_ref_start);
        for (int i = 0; i <= parameters.max_profile_array_size_int[j]; i++) {
            if(i < low_ref_start + exp_left_padding) {
                parameters.experimental_profile[j][i] = 0;
            } else if (i >= high_ref_end + exp_left_padding) {
                parameters.experimental_profile[j][i] = 1;
            } else {
                parameters.experimental_profile[j][i] = (parameters.experimental_profile[j][i] - low_ref) / (high_ref - low_ref);
            }
        }
    }
}

void
solver(double* control_parameters_mid, alglib::real_1d_array &residuals) // 
{
    //The 3 Concentration Arrays.
    double* D{new double[M]{}};
    double* ND{new double[M]{}};
    double* NS{new double[M]{}};

    double* D_out{new double[X * profile_file_line_count]{}};
    double* ND_out{new double[X * profile_file_line_count]{}};
    double* NS_out{new double[X * profile_file_line_count]{}};

    Eigen::MatrixXd m(X,X);
    Eigen::MatrixXd m_beads(X,X);
    Eigen::MatrixXd E(X,1);
    Eigen::MatrixXd E_NS(X,1);
    Eigen::MatrixXd E_ND(X,1);
    Eigen::MatrixXd solution(X,1);
    Eigen::MatrixXd solution_NS(X,1);
    Eigen::MatrixXd solution_ND(X,1);

    for(int i = 0; i < X; i++) {
        for (int j = 0; j < X; j++) {
            m(i,j) = (*matrixes.m)(i,j);
            m_beads(i,j) = (*matrixes.m_beads)(i,j);
            if (i == 0) {
                E(j) = (*matrixes.E)(j);
                E_NS(j) = (*matrixes.E_NS)(j);
                E_ND(j) = (*matrixes.E_ND)(j);
                solution(j) = (*matrixes.solution)(j);
                solution_NS(j) = (*matrixes.solution_NS)(j);
                solution_ND(j) = (*matrixes.solution_ND)(j);
            }
        }
    }
    
    double reaction_rate = 0.0d;
    double reverse_reaction_rate = 0.0d;
    for (int main_loop = 0; main_loop < profile_file_line_count; main_loop++) {
        //initialize the 3 Concentration Arrays. 
        for (int z = 0; z < Z; z++) {
            for (int x = 0; x < X; x++) {
                int i = x + X * z;

                if (x > (X / 2 - 1)) {
                    D[i]  = parameters.dye_conc[main_loop]; // molecules / m3
                    NS[i] = 0.0d;
                } else {
                    D[i]  = 0.0d;
                    NS[i] = parameters.bead_conc[main_loop];
                }
                ND[i] = 0.0d;
            }
        }

        for (int z = 1; z < Z; z++) {
            for (int i = 0; i < X; i++) {
                reaction_rate = control_parameters_mid[3] * dt * (D[i + (z - 1) * X] * NS[i + (z - 1) * X]); // - ND[i + (z - 1) * X] / Keq);
                reaction_rate = min(reaction_rate, NS[i + (z - 1) * X] / control_parameters_mid[2]);
                reaction_rate = min(D[i + (z - 1) * X], reaction_rate);
                reverse_reaction_rate = control_parameters_mid[4] * dt * ND[i + (z - 1) * X];
                reverse_reaction_rate = min(reverse_reaction_rate, ND[i + (z - 1) * X]);
                if (i == 0) {
                    E(i,0)      = (1.00d - r_dye)   * D[i + (z - 1) * X]    + r_dye   * D[i + (z - 1) * X + 1]                                                              - reaction_rate + reverse_reaction_rate;
                    E_NS(i,0)   = (1.00d - r_beads) * NS[i + (z - 1) * X]   + r_beads * NS[i + (z - 1) * X + 1]                                                             - reaction_rate / control_parameters_mid[2] + reverse_reaction_rate / control_parameters_mid[2];
                    E_ND(i,0)   = (1.00d - r_beads) * ND[i + (z - 1) * X]   + r_beads * ND[i + (z - 1) * X + 1]                                                             + reaction_rate - reverse_reaction_rate;
                    
                } else if (i == X - 1) {
                    E(i,0)      =                                             r_dye   * D[i + (z - 1) * X - 1]                  + (1.00d - r_dye)   * D[i + (z - 1) * X]    - reaction_rate + reverse_reaction_rate;
                    E_NS(i,0)   =                                             r_beads * NS[i + (z - 1) * X - 1]                 + (1.00d - r_beads) * NS[i + (z - 1) * X]   - reaction_rate / control_parameters_mid[2] + reverse_reaction_rate / control_parameters_mid[2];
                    E_ND(i,0)   =                                             r_beads * ND[i + (z - 1) * X - 1]                 + (1.00d - r_beads) * ND[i + (z - 1) * X]   + reaction_rate - reverse_reaction_rate;;
                } else {
                    E(i,0)      = r_dye   * D[i + (z - 1) * X - 1]          + (2.00d - 2.00d * r_dye)   * D[i + (z - 1) * X]    + r_dye   * D[i + (z - 1) * X + 1]          - reaction_rate + reverse_reaction_rate;
                    E_NS(i,0)   = r_beads * NS[i + (z - 1) * X - 1]         + (2.00d - 2.00d * r_beads) * NS[i + (z - 1) * X]   + r_beads * NS[i + (z - 1) * X + 1]         - reaction_rate / control_parameters_mid[2] + reverse_reaction_rate / control_parameters_mid[2];
                    E_ND(i,0)   = r_beads * ND[i + (z - 1) * X - 1]         + (2.00d - 2.00d * r_beads) * ND[i + (z - 1) * X]   + r_beads * ND[i + (z - 1) * X + 1]         + reaction_rate - reverse_reaction_rate;
                }
            }

            solution        = m * E;
            solution_NS     = m_beads * E_NS;
            solution_ND     = m_beads * E_ND;

            for (int i = 0; i < X; i++) {
                D[i + z * X]    = max(min(solution(i,0), parameters.dye_conc[main_loop]), 0.00d);
                D_out[i + main_loop * X]         = D[i + z * X];
                NS[i + z * X]   = max(min(solution_NS(i,0), parameters.bead_conc[main_loop]), 0.00d);
                NS_out[i + main_loop * X]        = NS[i + z * X];
                ND[i + z * X]   = max(solution_ND(i,0), 0.00d);
                ND_out[i + main_loop * X]        = ND[i + z * X];
            }
        }

        double split = 99.5d;
        int bottom_point;
        int top_point;
        num_profile_width = window_size - control_parameters_mid[0] - control_parameters_mid[1];    
        for (int i = 0; i < window_size; i++) {
            if (i < control_parameters_mid[0]) {
                parameters.numeric_model_profile[main_loop][i] = 0;
            } else if (i < window_size - control_parameters_mid[1]) {
                split = (double)(i - control_parameters_mid[0]) * ((double)X - 1.0d)/ (double)(num_profile_width - 1);
                bottom_point = static_cast<int>(floor(split));
                top_point = static_cast<int>(ceil(split));
                if (top_point == bottom_point) {
                    parameters.numeric_model_profile[main_loop][i] = D[bottom_point + (Z - 1) * X] + ND[bottom_point + (Z - 1) * X];
                } else {
                    parameters.numeric_model_profile[main_loop][i] = (D[bottom_point + (Z - 1) * X] + ND[bottom_point + (Z - 1) * X]) * ((double)top_point - split) + (D[top_point + (Z - 1) * X] + ND[top_point + (Z - 1) * X]) * (split - (double)bottom_point);
                }
            } else {
                parameters.numeric_model_profile[main_loop][i] = parameters.dye_conc[main_loop];
            }
            residuals[i + main_loop * window_size] = abs(parameters.numeric_model_profile[main_loop][i] / parameters.dye_conc[main_loop] - parameters.experimental_profile[main_loop][i]);
        }
    }
}

void
alglib_shift_solver(const alglib::real_1d_array &control_parameters, alglib::real_1d_array &residuals, void *ptr)
{
    double control_parameters_mid [5] = {control_parameters[0], control_parameters[1], parameters.p, parameters.kon, parameters.koff};
    solver(control_parameters_mid, residuals);
}

void
alglib_react_solver(const alglib::real_1d_array &react_control_parameters, alglib::real_1d_array &residuals, void *ptr)
{
    double control_parameters_mid [5] = {(double)num_left_padding, (double)num_right_padding, react_control_parameters[0], react_control_parameters[1], react_control_parameters[2]};
    solver(control_parameters_mid, residuals);
}

void
alglib_solver(const alglib::real_1d_array &control_parameters, alglib::real_1d_array &residuals, void *ptr)
{
    double control_parameters_mid [5] = {control_parameters[0], control_parameters[1], control_parameters[2], control_parameters[3], control_parameters[4]};
    solver(control_parameters_mid, residuals);
}

int
main()
{
    parameters.wt_percent                   = new double[profile_file_line_count];   //= 0.1d; // wt%
    parameters.dye_conc_mgml                = new double[profile_file_line_count];   //= 0.00336d; // mg/ml FITC
    parameters.bead_conc                    = new double[profile_file_line_count];   //= wt_percent / 100.0d / 1.05d / ( 4.0d / 3.0d * M_PI * powf(diameter / 2.0d, 3.0d)); // beads/m3 
    parameters.dye_conc                     = new double[profile_file_line_count];   //= dye_conc_mgml * 1000.0d / 332.326d * 6.022e+23;     // molecules FITC / m3
    parameters.beginning_of_channel         = new int[profile_file_line_count];     //= 56;
    parameters.end_of_channel               = new int[profile_file_line_count];     //= 196;
    parameters.profile_array_size_int       = new int[profile_file_line_count];     //= 140;
    parameters.max_profile_array_size_int   = new int[profile_file_line_count];

    //The 3 Concentration Arrays.
    double* D{new double[M]{}};
    double* ND{new double[M]{}};
    double* NS{new double[M]{}};

    double* D_out{new double[X * profile_file_line_count]{}};
    double* ND_out{new double[X * profile_file_line_count]{}};
    double* NS_out{new double[X * profile_file_line_count]{}};

    Eigen::MatrixXd m           = Eigen::MatrixXd::Zero(X, X);
    matrixes.m                  = &m;
    Eigen::MatrixXd m_beads     = Eigen::MatrixXd::Zero(X, X);
    matrixes.m_beads            = &m_beads;
    Eigen::MatrixXd E           = Eigen::MatrixXd::Zero(1, X);
    matrixes.E                  = &E;
    Eigen::MatrixXd E_NS        = Eigen::MatrixXd::Zero(1, X);
    matrixes.E_NS               = &E_NS;
    Eigen::MatrixXd E_ND        = Eigen::MatrixXd::Zero(1, X);
    matrixes.E_ND               = &E_ND;
    Eigen::MatrixXd solution    = Eigen::MatrixXd::Zero(1, X);
    matrixes.solution           = &solution;
    Eigen::MatrixXd solution_NS = Eigen::MatrixXd::Zero(1, X);
    matrixes.solution_NS        = &solution_NS;
    Eigen::MatrixXd solution_ND = Eigen::MatrixXd::Zero(1, X);
    matrixes.solution_ND        = &solution_ND;

    for (int i = 0; i < X; i++) {
        for (int j = 0; j < X; j++) {
            if (j == i - 1) {
                m(i , j) = -r_dye;
                m_beads(i , j) = -r_beads;
            } else if (j == i + 1) {
                m(i , j) = -r_dye;
                m_beads(i , j) = -r_beads;
            } else if ((i != 0 && i != X - 1) && (i == j)) {
                m(i , j) = 2.0d + 2.0d * r_dye;
                m_beads(i , j) = 2.0d + 2.0d * r_beads;
            } else if ((i == 0 || i == X - 1) && (i == j)) {
                m(i , j) = 1.0d + r_dye;
                m_beads(i , j) = 1.0d + r_beads;
            } else {
                m(i , j) = 0.0d;
                m_beads(i , j) = 0.0d;
            }
        }
    }
    m = m.inverse();
    m_beads = m_beads.inverse();

    for (int i = 0; i < X * profile_file_line_count; i++) {
        D_out[i] = 0;
        ND_out[i] = 0;
        NS_out[i] = 0;
    }

    std::cout << "Priming " << output_file_name << ": ";
    prime_excel_output(output_file_name);
    std::cout << "Done" << endl;

    std::cout << "Reading " << input_file_name << ": ";
    read_profile_from_text(input_file_name);
    std::cout << "Done" << endl;

    std::cout << "Normalizing " << profile_file_line_count << " profiles: ";
    normalize_profile();
    std::cout << "Done" << endl;

    try
    {
        double epsx = 0.0000000001;
        alglib::ae_int_t maxits = 0;
        //alglib::real_1d_array residuals[profile_file_line_count * window_size];

        double init_shift_control_parameters[] = {(double)num_left_padding, (double)num_right_padding};
        alglib::real_1d_array shift_control_parameters;
        shift_control_parameters.setcontent(2, init_shift_control_parameters);
        /*
        alglib::real_1d_array shift_s = "[1.0, 1.0]";
        alglib::real_1d_array shift_bndl = "[0, 0]";
        alglib::real_1d_array shift_bndu = "[60, 60]";
        alglib::minlmstate shift_state;
        alglib::minlmreport shift_rep;
        alglib::minlmcreatev(2, window_size * profile_file_line_count, shift_control_parameters, 0.0001, shift_state);
        alglib::minlmsetbc(shift_state, shift_bndl, shift_bndu);
        alglib::minlmsetcond(shift_state, epsx, maxits);
        alglib::minlmsetscale(shift_state, shift_s);
        alglib::minlmsetnonmonotonicsteps(shift_state, 2);
        std::cout << "shift minlmoptimize: ";
        alglib::minlmoptimize(shift_state, alglib_shift_solver);   // Optimize
        std::cout << "Done" << endl;
        alglib::minlmresults(shift_state, shift_control_parameters, shift_rep);
        num_left_padding = shift_control_parameters[0];
        num_right_padding = shift_control_parameters[1];
        printf("%s\n", shift_control_parameters.tostring(4).c_str());
        alglib::minlmrestartfrom(shift_state, shift_control_parameters);
*/
        double init_react_control_parameters[] = {parameters.p, parameters.kon, parameters.koff};
        alglib::real_1d_array react_control_parameters;
        react_control_parameters.setcontent(3, init_react_control_parameters);
        /*
        alglib::real_1d_array react_s = "[0.1,1.0e-19,1.0e-19]";
        alglib::real_1d_array react_bndl = "[10,1.0e-24,1.0e-24]";
        alglib::real_1d_array react_bndu = "[200,1.0e-1,1.0e-16]";
        alglib::minlmstate react_state;
        alglib::minlmreport react_rep;
        alglib::minlmcreatev(3, window_size * profile_file_line_count, react_control_parameters, 0.000001, react_state);
        alglib::minlmsetbc(react_state, react_bndl, react_bndu);
        alglib::minlmsetcond(react_state, epsx, maxits);
        alglib::minlmsetscale(react_state, react_s);
        alglib::minlmsetnonmonotonicsteps(react_state, 2);
        std::cout << "react minlmoptimize: ";
        alglib::minlmoptimize(react_state, alglib_react_solver);   // Optimize
        std::cout << "Done" << endl;
        alglib::minlmresults(react_state, react_control_parameters, react_rep);
        printf("%s\n", react_control_parameters.tostring(4).c_str());
        */
        double init_control_parameters[] = {shift_control_parameters[0], shift_control_parameters[1], react_control_parameters[2], react_control_parameters[3]};
        alglib::real_1d_array control_parameters;
        control_parameters.setcontent(5, init_control_parameters);
        alglib::real_1d_array s = "[1.0, 1.0, 0.1, 1.0e-18, 1.0e-18]";
        alglib::real_1d_array bndl = "[0, 0, 10, 1.0e-22, 0]";
        alglib::real_1d_array bndu = "[60, 60, 1000, 1.0e-17, 1.0e-17]";
        alglib::minlmstate state;
        alglib::minlmreport rep;
        alglib::minlmcreatev(5, window_size * profile_file_line_count, control_parameters, 0.000001, state);
        alglib::minlmsetbc(state, bndl, bndu);
        alglib::minlmsetcond(state, epsx, maxits);
        alglib::minlmsetscale(state, s);
        alglib::minlmsetnonmonotonicsteps(state, 2);
        std::cout << "minlmoptimize: ";
        alglib::minlmoptimize(state, alglib_solver);   // Optimize
        std::cout << "Done" << endl;
        alglib::minlmresults(state, control_parameters, rep);
        printf("%s\n", control_parameters.tostring(4).c_str());

        double control_parameters_out[5] = {control_parameters[0], control_parameters[1], control_parameters[2], control_parameters[3], control_parameters[4]};   
        save_excel_output(D_out, ND_out, NS_out, control_parameters_out);
        
        //printf("%s\n", control_parameters.tostring(2).c_str()); // EXPECTED: [-3,+3]
    }
    catch(alglib::ap_error alglib_exception)
    {
        printf("ALGLIB exception with message '%s'\n", alglib_exception.msg.c_str());
        return 1;
    }
    
    std::cout << "Program finished " << endl;
    return 0;
}