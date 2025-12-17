#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include <eigen-3.4.0/Eigen/Dense>
// D: Dye
// NS: Nanoparticle Site
// ND: Nanoparticle Bound Dye

using namespace std;

void prime_excel_output(string file_name);
void save_excel_output(string array_name, double* double_array);
void solver(); 
void set_matrixes(double new_flow, double new_diameter);
void save_variable_excel_output(string array_name, double* double_array, int array_length);

const int Z                 = 100; // divisions along the z(time) axis
const int X                 = 500; // divisions along the x axis
const int M                 = X * Z;

// Device Dimenssions
const double W = 5e-4;  //meters: 500 um
const double H = 4e-5;  //meters: 40 um
const double L = 0.025; //meters: 2.5 cm

// Operating conditions/settings
double flow            = 2.0d * 5e-9d / 60;                        //m3/s: 2*5 ulmin
double diameter        = 41e-9d; // meters
const int window_size = 180;

// Physical Constants
const double visc            = 0.0010016d; // Dynamic viscosity of water at 20C in Pa.s
const double difusion_dye    = 4.9e-10d; // m2/s From 4.9 × 10−6 cm2 s−1 The diffusion coefficient of fluorescein in water at 21.5°C, as calculated from the Wilke-Chang correlation
double difusion_beads  = 1.380649e-23d * (273.15d + 25.0d) / (3.0d * M_PI * visc * diameter); // m2/s From kB*T/(3*pi*visc*d) kB=1.380649×10−23 J⋅K−1

// Derived values for Crank-Nicolson implicit method
double restime       = W * H * L / flow; // seconds
double dt            = restime / Z; // seconds
double T_dye         = difusion_dye * restime / (W * W);
double T_beads       = difusion_beads * restime / (W * W);
double dT_dye        = difusion_dye * dt / (W * W); // T = Dt/l^2
double dT_beads      = difusion_beads * dt / (W * W);
const double dx            = W / X; // m / subdivision
const double dX            = 1.0d / X; // X = x/l
double r_dye         = dT_dye / dX / dX; // r = dT/(dX)^2
double r_beads       = dT_beads / dX / dX; // r = dT/(dX)^2

struct parameters_struct {
    double p               = 50.3152d; // FITC molecules / PS bead
    double kon             = 2.4526E-19; // 0.10d;
    double koff            = 11.80657d;//1.987E-16d * (W / X);
    double* scale;
    double wt_percent;                      //= 0.1d; // wt%
    double dye_conc_mgml;                   //= 0.00336d; // mg/ml FITC
    double bead_conc;                       //= wt_percent / 100.0d / 1.05d / ( 4.0d / 3.0d * M_PI * powf(diameter / 2.0d, 3.0d)); // beads/m3 
    double dye_conc;                        //= dye_conc_mgml * 1000.0d / 332.326d * 6.022e+23;     // molecules FITC / m3
    double experimental_profile[window_size];    // 15_111 16_140
    double numeric_model_profile[window_size];
    double* D;
    double* NS;
    double* ND;
    double* D_out;
    double* NS_out;
    double* ND_out;
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
set_matrixes(double new_flow, double new_diameter)
{
    flow            = new_flow;                        //m3/s: 2*5 ulmin
    diameter        = new_diameter; // meters
    difusion_beads  = 1.380649e-23d * (273.15d + 25.0d) / (3.0d * M_PI * visc * diameter); // m2/s From kB*T/(3*pi*visc*d) kB=1.380649×10−23 J⋅K−1
    restime         = W * H * L / flow; // seconds
    dt              = restime / Z; // seconds
    T_dye           = difusion_dye * restime / (W * W);
    T_beads         = difusion_beads * restime / (W * W);
    dT_dye          = difusion_dye * dt / (W * W); // T = Dt/l^2
    dT_beads        = difusion_beads * dt / (W * W);
    r_dye           = dT_dye / dX / dX; // r = dT/(dX)^2
    r_beads         = dT_beads / dX / dX; // r = dT/(dX)^2

    Eigen::MatrixXd m(X,X);
    for(int i = 0; i < X; i++) {
        for (int j = 0; j < X; j++) {
            m(j,i) = (*matrixes.m)(j,i);
        }
    }
    Eigen::MatrixXd m_beads(X,X);
    for(int i = 0; i < X; i++) {
        for (int j = 0; j < X; j++) {
            m_beads(j,i) = (*matrixes.m_beads)(j,i);
        }
    }

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

    for(int i = 0; i < X; i++) {
        for (int j = 0; j < X; j++) {
            (*matrixes.m)(j,i) = m(j,i);
            (*matrixes.m_beads)(j,i) = m_beads(j,i);
        }
    }
}

void
prime_excel_output(string file_name)
{
    ofstream fout;
    fout.open(file_name, std::ofstream::out | std::ofstream::trunc);
    fout << "res_time  " << "bind_ratio(p)  " << "forward_reaction_rate  "  << "reverse_reaction_rate  "<< "dye_conc.   " << "bead_conc.  " << "bead_diameter   "; 
    for (int x = 0; x < X; x++) {
        fout << x << "  ";
    }
    fout << endl;
    fout << "sec  " << "molc/bead  " << "forward_reaction_rate  " << "reverse_reaction_rate  " << "umol/l   " << "10E+18 bead/m3  " << "um    ";
    fout << endl;
    fout.close(); 
}

void
save_excel_output(string array_name, double* double_array)
{
    ofstream fout;
    fout.open(array_name, std::ofstream::out | std::ofstream::app);
    fout << (Z - 1) * dt << "  " << parameters.p << "  " << parameters.kon << "  " << parameters.koff << "  " << parameters.dye_conc / (6.02214076e23) * 1000 << "   " <<  parameters.bead_conc / 1.0e18 << "  " << diameter / 1e-9 << "  ";
    for (int x = 0; x < X; x++) {
        fout << double_array[x] << "    ";
    }
    fout << endl;
    fout.close(); 
}

void
save_variable_excel_output(string array_name, double* double_array, int array_length)
{
    ofstream fout;
    fout.open(array_name, std::ofstream::out | std::ofstream::app);
    fout << (Z - 1) * dt << "  " << parameters.p << "  " << parameters.kon << "  " << parameters.koff << "  " << parameters.dye_conc / (6.02214076e23) * 1000 << "   " <<  parameters.bead_conc / 1.0e18 << "  " << diameter / 1e-9 << "  ";
    for (int x = 0; x < array_length; x++) {
        fout << double_array[x] << "    ";
    }
    fout << endl;
    fout.close(); 
}

void
solver() // 
{
    //The 3 Concentration Arrays.
    double* D = parameters.D;
    double* ND = parameters.ND;
    double* NS = parameters.NS; //{new double[M]{}}

    double* D_out = parameters.D_out;//{new double[X * profile_file_line_count]{}};
    double* ND_out = parameters.ND_out;//{new double[X * profile_file_line_count]{}};
    double* NS_out = parameters.NS_out; //{new double[X * profile_file_line_count]{}};

    Eigen::MatrixXd m(X,X);
    for(int i = 0; i < X; i++) {
        for (int j = 0; j < X; j++) {
            m(j,i) = (*matrixes.m)(j,i);
        }
    }
    Eigen::MatrixXd m_beads(X,X);
    for(int i = 0; i < X; i++) {
        for (int j = 0; j < X; j++) {
            m_beads(j,i) = (*matrixes.m_beads)(j,i);
        }
    }
    Eigen::MatrixXd E(X,1);
    for (int j = 0; j < X; j++) {
        E(j) = (*matrixes.E)(j);
    }
    Eigen::MatrixXd E_NS(X,1);
    for (int j = 0; j < X; j++) {
        E_NS(j) = (*matrixes.E_NS)(j);
    }      
    Eigen::MatrixXd E_ND(X,1);
    for (int j = 0; j < X; j++) {
        E_ND(j) = (*matrixes.E_ND)(j);
    }       
    Eigen::MatrixXd solution(X,1);
    for (int j = 0; j < X; j++) {
        solution(j) = (*matrixes.solution)(j);
    }     
    Eigen::MatrixXd solution_NS(X,1);
    for (int j = 0; j < X; j++) {
        solution_NS(j) = (*matrixes.solution_NS)(j);
    }    
    Eigen::MatrixXd solution_ND(X,1);
    for (int j = 0; j < X; j++) {
        solution_ND(j) = (*matrixes.solution_ND)(j);
    }       
    double reaction_rate = 0.0d;
    
    for (int z = 0; z < Z; z++) {
        for (int x = 0; x < X; x++) {
            int i = x + X * z;

            if (x > (X / 2 - 1)) {
                D[i]  = parameters.dye_conc; // molecules / m3
                NS[i] = 0.0d;
            } else {
                D[i]  = 0.0d;
                NS[i] = parameters.bead_conc;
            }
            ND[i] = 0.0d;
        }
    }

    for (int z = 1; z < Z; z++) {
        for (int i = 0; i < X; i++) {
            reaction_rate = parameters.kon * dt * (D[i + (z - 1) * X] * NS[i + (z - 1) * X]) - parameters.koff * dt * ND[i + (z - 1) * X]; // - ND[i + (z - 1) * X] / Keq);
            reaction_rate = min(reaction_rate, NS[i + (z - 1) * X] * parameters.p);
            reaction_rate = min(D[i + (z - 1) * X], reaction_rate);
            reaction_rate = max(-ND[i + (z - 1) * X], reaction_rate);
            if (i == 0) {
                E(i,0)      = (1.00d - r_dye)   * D[i + (z - 1) * X]    + r_dye   * D[i + (z - 1) * X + 1]                                                              - reaction_rate;
                E_NS(i,0)   = (1.00d - r_beads) * NS[i + (z - 1) * X]   + r_beads * NS[i + (z - 1) * X + 1]                                                             - reaction_rate / parameters.p;
                E_ND(i,0)   = (1.00d - r_beads) * ND[i + (z - 1) * X]   + r_beads * ND[i + (z - 1) * X + 1]                                                             + reaction_rate;
                
            } else if (i == X - 1) {
                E(i,0)      =                                             r_dye   * D[i + (z - 1) * X - 1]                  + (1.00d - r_dye)   * D[i + (z - 1) * X]    - reaction_rate;
                E_NS(i,0)   =                                             r_beads * NS[i + (z - 1) * X - 1]                 + (1.00d - r_beads) * NS[i + (z - 1) * X]   - reaction_rate / parameters.p;
                E_ND(i,0)   =                                             r_beads * ND[i + (z - 1) * X - 1]                 + (1.00d - r_beads) * ND[i + (z - 1) * X]   + reaction_rate;
            } else {
                E(i,0)      = r_dye   * D[i + (z - 1) * X - 1]          + (2.00d - 2.00d * r_dye)   * D[i + (z - 1) * X]    + r_dye   * D[i + (z - 1) * X + 1]          - reaction_rate;
                E_NS(i,0)   = r_beads * NS[i + (z - 1) * X - 1]         + (2.00d - 2.00d * r_beads) * NS[i + (z - 1) * X]   + r_beads * NS[i + (z - 1) * X + 1]         - reaction_rate / parameters.p;
                E_ND(i,0)   = r_beads * ND[i + (z - 1) * X - 1]         + (2.00d - 2.00d * r_beads) * ND[i + (z - 1) * X]   + r_beads * ND[i + (z - 1) * X + 1]         + reaction_rate;
            }
        }

        solution        = m * E;
        solution_NS     = m_beads * E_NS;
        solution_ND     = m_beads * E_ND;

        for (int i = 0; i < X; i++) {
            D[i + z * X]    = max(min(solution(i,0), parameters.dye_conc), 0.00d);
            D_out[i]         = D[i + z * X];
            NS[i + z * X]   = max(min(solution_NS(i,0), parameters.bead_conc), 0.00d);
            NS_out[i]        = NS[i + z * X];
            ND[i + z * X]   = max(solution_ND(i,0), 0.00d);
            ND_out[i]        = ND[i + z * X];
        }
    }

    double split = 99.5d;
    int bottom_point;
    int top_point;  
    for (int i = 0; i < window_size; i++) {
        if (i < 0.0d) {
            parameters.numeric_model_profile[i] = 0;
        } else if (i < window_size) {
            split = (double)(i) * ((double)X - 1.0d) / (double)(window_size - 1.0d);
            bottom_point = static_cast<int>(floor(split));
            top_point = static_cast<int>(ceil(split));
            if (top_point == bottom_point || bottom_point == (X - 1)) {
                parameters.numeric_model_profile[i] = D[bottom_point + (Z - 1) * X] + ND[bottom_point + (Z - 1) * X];
            } else {
                parameters.numeric_model_profile[i] = (D[bottom_point + (Z - 1) * X] + ND[bottom_point + (Z - 1) * X]) * ((double)top_point - split) + (D[top_point + (Z - 1) * X] + ND[top_point + (Z - 1) * X]) * (split - (double)bottom_point);
            }
        } else {
            parameters.numeric_model_profile[i] = parameters.dye_conc;
        }
    }
}

int
main()
{
    //The 3 Concentration Arrays.
    double* D{new double[M]{}};
    double* ND{new double[M]{}};
    double* NS{new double[M]{}};
    parameters.D = D;
    parameters.ND = ND;
    parameters.NS = NS;
    
    double* D_out{new double[X]{}};
    double* ND_out{new double[X]{}};
    double* NS_out{new double[X]{}};
    double* dye_total{new double[X]{}};
    double* bead_total{new double[X]{}};
    double* bead_bound{new double[X]{}};
    double* dye_total_norm{new double[X]{}};
    double* derivative{new double[X]{}};
    double* d_a{new double[16]{}};
    double min_a = 0;
    double max_d = 0;
    parameters.D_out = D_out;
    parameters.ND_out = ND_out;
    parameters.NS_out = NS_out;

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

    set_matrixes(flow, diameter);

    for (int i = 0; i < X; i++) {
        D_out[i] = 0;
        ND_out[i] = 0;
        NS_out[i] = 0;
    }

    std::cout << "Priming arrays: ";
    /*
    prime_excel_output("dye_total_array.txt");
    prime_excel_output("dye_free_array.txt");
    prime_excel_output("dye_bound_array.txt");
    prime_excel_output("bead_total_array.txt");
    prime_excel_output("bead_free_array.txt");
    prime_excel_output("bead_bound_array.txt");
    prime_excel_output("normalized_dye_total.txt");
    prime_excel_output("derivative.txt");
    */
    prime_excel_output("DA.txt");
    std::cout << "Done" << endl;
    
    for (int scan = 0; scan < 8; scan++) {
        std::cout << "Running Model: ";
        switch(scan) {
            case 0:
                parameters.p = 109.3;
                parameters.kon = 4.07926e-19d;
                diameter = 41e-9d;
                parameters.koff = 0.0;
                parameters.wt_percent = 0.08;
            break;
            case 1:
                parameters.p = 109.3 * 0.5;
                parameters.kon = 4.07926e-19d;
                diameter = 41e-9d;
                parameters.koff = 0.0;
                parameters.wt_percent = 0.08;
            break;
            case 2:
                parameters.p = 109.3 * 2.0;
                parameters.kon = 4.07926e-19d;
                diameter = 41e-9d;
                parameters.koff = 0.0;
                parameters.wt_percent = 0.08;
            break;
            case 3:
                parameters.p = 109.3;
                parameters.kon = 0.5 * 4.07926e-19d;
                diameter = 41e-9d;
                parameters.koff = 0.0;
                parameters.wt_percent = 0.08;
            break;
            case 4:
                parameters.p = 109.3;
                parameters.kon = 2.0 * 4.07926e-19d;
                diameter = 41e-9d;
                parameters.koff = 0.0;
                parameters.wt_percent = 0.08;
            break;
            case 5:
                parameters.p = 109.3;
                parameters.kon = 4.07926e-19d;
                diameter = 0.5 * 41e-9d;
                parameters.koff = 0.0;
                parameters.wt_percent = 0.08;
            break;
            case 6:
                parameters.p = 109.3;
                parameters.kon = 4.07926e-19d;
                diameter = 2.0 * 41e-9d;
                parameters.koff = 0.0;
                parameters.wt_percent = 0.08;
            break;
            case 7:
                parameters.p = 109.3;
                parameters.kon = 4.07926e-19d;
                diameter = 41e-9d;
                parameters.koff = 10.0;
                parameters.wt_percent = 0.08;
            break;
        }
        
        flow;
        for(int wt_scan = 0; wt_scan < 16; wt_scan++) {
            parameters.wt_percent = 0.01 * wt_scan;
            parameters.bead_conc = parameters.wt_percent / 100.0d / 1.05d / ( 4.0d / 3.0d * M_PI * powf(diameter / 2.0d, 3.0d)); // beads/m3 
            parameters.dye_conc_mgml = 0.0033d;
            parameters.dye_conc = parameters.dye_conc_mgml * 1000.0d / 332.326d * 6.022e+23;     // molecules FITC / m3

            set_matrixes(flow, diameter);
            solver();

            for (int i = 0; i < X; i++) {
                dye_total[i] = D_out[i] + ND_out[i];
                bead_total[i] = NS_out[i] + ND_out[i] / parameters.p;
                bead_bound[i] = ND_out[i] / parameters.p;
                dye_total_norm[i] = dye_total[i] / parameters.dye_conc;
            }
            min_a = 0;
            max_d = 0;
            for (int i =0; i < X; i++) {
                if(i == 0 || i == (X-1)) {
                    derivative[i] = 0;
                } else {
                    derivative[i] = (dye_total_norm[i+1] - dye_total_norm[i-1]) / 2; 
                }
                min_a = min(min_a, derivative[i]);
                max_d = max(max_d, derivative[i]);
            }
            d_a[wt_scan] = max_d - min_a;
        }
        save_variable_excel_output("DA.txt", d_a, 16);
        /*
        save_excel_output("dye_total_array.txt", dye_total);
        save_excel_output("dye_free_array.txt", D_out);
        save_excel_output("dye_bound_array.txt", ND_out);
        save_excel_output("bead_total_array.txt", bead_total);
        save_excel_output("bead_free_array.txt", NS_out);
        save_excel_output("bead_bound_array.txt", bead_bound);
        save_excel_output("normalized_dye_total.txt", dye_total_norm);
        save_excel_output("derivative.txt", derivative);
        */

        std::cout << "Done" << endl;
    }

    std::cout << "Program finished " << endl;
    return 0;
}