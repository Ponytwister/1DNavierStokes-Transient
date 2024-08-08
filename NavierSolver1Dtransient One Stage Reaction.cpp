#include <cmath>
#include <fstream>
#include <iostream>
#include <eigen-3.4.0/Eigen/Dense>

// D: Dye
// NS: Nanoparticle Site
// ND: Nanoparticle Bound Dye

using namespace std;
void SetInletCondit(float* D, float* ND, float* NS, int z = 0);

void SaveExcelOutput(int z, float* D, float* ND, float* NS, float p, float kon);

const int Z                 = 100; // divisions along the z(time) axis
const int X                 = 400; // divisions along the x axis
const int M                 = X * Z;

// Device Dimensions
const float W = 5e-4;  //meters: 500 um
const float H = 4e-5;  //meters: 40 um
const float L = 0.025; //meters: 2.5 cm

// Operating conditions/settings
const float flow            = 2.0f * 5e-9f / 60;                        //m3/s: 2*5 ulmin
const float diameter        = 41e-9f; // meters
const float wt_percent      = 0.1f; // wt%
const float dye_conc_mgml   = 0.00336f; // mg/ml FITC
const float bead_conc             = wt_percent / 100.0f / 1.05f / ( 4.0f / 3.0f * M_PI * powf(diameter / 2.0f, 3.0f)); // beads/m3 
const float dye_conc              = dye_conc_mgml * 1000.0f / 332.326f * 6.022e+23;     // molecules FITC / m3 

// Physical Constants
const float visc            = 0.0010016f; // Dynamic viscosity of water at 20C in Pa.s
const float difusion_dye    = 4.9e-10f; // m2/s From 4.9 × 10−6 cm2 s−1 The diffusion coefficient of fluorescein in water at 21.5°C, as calculated from the Wilke-Chang correlation
const float difusion_beads  = 1.380649e-23f * (273.15f + 25.0f) / (3.0f * M_PI * visc * diameter); // m2/s From kB*T/(3*pi*visc*d) kB=1.380649×10−23 J⋅K−1

// Derived values for Crank-Nicolson implicit method
const float restime       = W * H * L / flow; // seconds
const float dt            = restime / Z; // seconds
const float T_dye         = difusion_dye * restime / (W * W);
const float T_beads       = difusion_beads * restime / (W * W);
const float dT_dye        = difusion_dye * dt / (W * W); // T = Dt/l^2
const float dT_beads      = difusion_beads * dt / (W * W);
const float dx            = W / X; // m / subdivision
const float dX            = 1.0f / X; // X = x/l
const float r_dye         = dT_dye / dX / dX; // r = dT/(dX)^2
const float r_beads       = dT_beads / dX / dX; // r = dT/(dX)^2

int
main()
{
    std::cout << "restime = "   << restime  << endl;
    std::cout << "dt = "        << dt       << endl;
    std::cout << "T_dye = "     << T_dye    << endl;
    std::cout << "T_beads = "   << T_beads  << endl;
    std::cout << "dT_dye = "    << dT_dye   << endl;
    std::cout << "dT_beads = "  << dT_beads << endl;
    std::cout << "dx = "        << dx       << endl;
    std::cout << "dX = "        << dX       << endl;
    std::cout << "r_dye = "     << r_dye    << endl;
    std::cout << "r_beads = "   << r_beads  << endl;

    //The 3 Concentration Arrays.
    float* D{new float[M]{}};
    float* ND{new float[M]{}};
    float* NS{new float[M]{}};

    //initialize the 3 Concentration Arrays.
    for (int z = 0; z < Z; z++) {
        SetInletCondit(D, ND, NS, z);
    }

    Eigen::MatrixXf m(X,X);
    m.setZero();
    Eigen::MatrixXf m_beads(X,X);
    m_beads.setZero();
    Eigen::MatrixXf E(X,1);
    E.setZero();
    Eigen::MatrixXf E_NS(X,1);
    E_NS.setZero();
    Eigen::MatrixXf E_ND(X,1);
    E_ND.setZero();
    Eigen::MatrixXf Solution(X,1);
    Solution.setZero();
    Eigen::MatrixXf Solution_NS(X,1);
    Solution_NS.setZero();
    Eigen::MatrixXf Solution_ND(X,1);
    Solution_ND.setZero();


    for (int i = 0; i < X; i++) {
        for (int j = 0; j < X; j++) {
            if (j == i - 1) {
                m(i , j) = -r_dye;
                m_beads(i , j) = -r_beads;
            } else if (j == i + 1) {
                m(i , j) = -r_dye;
                m_beads(i , j) = -r_beads;
            } else if ((i != 0 && i != X - 1) && (i == j)) {
                m(i , j) = 2.0f + 2.0f * r_dye;
                m_beads(i , j) = 2.0f + 2.0f * r_beads;
            } else if ((i == 0 || i == X - 1) && (i == j)) {
                m(i , j) = 1.0f + r_dye;
                m_beads(i , j) = 1.0f + r_beads;
            } else {
                m(i , j) = 0.0f;
                m_beads(i , j) = 0.0f;
            }
        }
    }
    m = m.inverse();
    m_beads = m_beads.inverse();

    float P1              = 100.0f; // FITC molecules / PS bead
    //const float Keq       = 1422.0f;
    float kon1            = 1.0E-10f * (W / X); // 0.10f;
    float reaction_rate_1 = 0.0f;

        for (int z = 1; z < Z; z++) {
            for (int i = 0; i < X; i++) {
                reaction_rate_1 = kon1 * dt * (D[i + (z - 1) * X] * NS[i + (z - 1) * X]); // - ND[i + (z - 1) * X] / Keq);
                reaction_rate_1 = min(reaction_rate_1, NS[i + (z - 1) * X] * P1);
                reaction_rate_1 = min(D[i + (z - 1) * X], reaction_rate_1);
                if (i == 0) {
                    E(i,0)      = (1.00f - r_dye)   * D[i + (z - 1) * X]    + r_dye   * D[i + (z - 1) * X + 1]                                                              - reaction_rate_1;
                    E_NS(i,0)   = (1.00f - r_beads) * NS[i + (z - 1) * X]   + r_beads * NS[i + (z - 1) * X + 1]                                                             - reaction_rate_1 * P1;
                    E_ND(i,0)   = (1.00f - r_beads) * ND[i + (z - 1) * X]   + r_beads * ND[i + (z - 1) * X + 1]                                                             + reaction_rate_1;
                    
                } else if (i == X - 1) {
                    E(i,0)      =                                           r_dye   * D[i + (z - 1) * X - 1]                    + (1.00f - r_dye)   * D[i + (z - 1) * X]    - reaction_rate_1;
                    E_NS(i,0)   =                                           r_beads * NS[i + (z - 1) * X - 1]                   + (1.00f - r_beads) * NS[i + (z - 1) * X]   - reaction_rate_1 * P1;
                    E_ND(i,0)   =                                           r_beads * ND[i + (z - 1) * X - 1]                   + (1.00f - r_beads) * ND[i + (z - 1) * X]   + reaction_rate_1;
                } else {
                    E(i,0)      = r_dye   * D[i + (z - 1) * X - 1]          + (2.00f - 2.00f * r_dye)   * D[i + (z - 1) * X]    + r_dye   * D[i + (z - 1) * X + 1]          - reaction_rate_1;
                    E_NS(i,0)   = r_beads * NS[i + (z - 1) * X - 1]         + (2.00f - 2.00f * r_beads) * NS[i + (z - 1) * X]   + r_beads * NS[i + (z - 1) * X + 1]         - reaction_rate_1 / P1;
                    E_ND(i,0)   = r_beads * ND[i + (z - 1) * X - 1]         + (2.00f - 2.00f * r_beads) * ND[i + (z - 1) * X]   + r_beads * ND[i + (z - 1) * X + 1]         + reaction_rate_1;
                }
            }

            Solution        = m*E;
            Solution_NS     = m_beads*E_NS;
            Solution_ND     = m_beads*E_ND;


            for (int i = 0; i < X; i++) {
                D[i + z * X]    = max(min(Solution(i,0), dye_conc), 0.00f);
                NS[i + z * X]   = max(min(Solution_NS(i,0), bead_conc), 0.00f);
                ND[i + z * X]   = max(Solution_ND(i,0), 0.00f);
            }
            SaveExcelOutput(z-1, D, ND, NS, P1, kon1);
        }
    SaveExcelOutput(Z-1, D, ND, NS, P1, kon1);
    
    std::cout << "Program finished " << endl;
    return 0;
}

void
SetInletCondit(float* D, float* ND, float* NS, int z)
{
    int i = 0;
    for (int x = 0; x < X; x++) {
        i = x + X * z;

        if (x > (X / 2 - 1)) {
            D[i]  = dye_conc; // molecules / m3
            NS[i] = 0.0f;
        } else {
            D[i]  = 0.0f;
            NS[i] = bead_conc;
        }
        ND[i] = 0.0f;
    }
}


void
SaveExcelOutput(int z, float* D, float* ND, float* NS, float p, float kon)
{
    float Sum;
    float* Point;
    ofstream fout;
    fout.open("ExcelProfiles.txt", std::ofstream::out | std::ofstream::app);
    for (int j = 0; j < 3; j++) {
        fout << z * dt << "  " << p << "  " << kon << "  " << dye_conc << "   " <<  bead_conc << "  ";
        if (j == 0) {
            Point = D;
            fout << "D  ";
        }
        if (j == 1) {
            Point = ND;
            fout << "ND  ";
        }
        if (j == 2) {
            Point = NS;
            fout << "NS  ";
        }
        for (int x = 0; x < X; x++) {
            Sum = Point[x + z * X];
            fout << Sum << "    ";
        }
        fout << endl;
    }
    fout.close();
}