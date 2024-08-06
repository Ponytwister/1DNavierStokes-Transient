#include <cmath>
#include <fstream>
#include <iostream>
#include <eigen-3.4.0/Eigen/Dense>

// D: Dye
// NS: Nanoparticle Site
// ND: Nanoparticle Bound Dye

using namespace std;
void SetInletCondit(float* D, float* DD, float* ND, float* NS, float* NDD, int z = 0);
void SaveViewable(float* V, int X, int Y, int Z, string Filename);
void SaveExcelOutput(int z, float* D, float* DD, float* ND, float* NDD, float* NS);

const int Z                 = 100;
const int X                 = 400;
const int M                 = X * Z;

// Device Dimensions
const float W = 5e-4;  //meters: 500 um
const float H = 4e-5;  //meters: 40 um
const float L = 0.025; //meters: 2.5 cm

// Operating conditions/settings
const float Flow            = 2.0f * 5e-9f / 60;                        //m3/s: 2*5 ulmin
const float diameter        = 23e-9f;
float BeadConc              = 0.001f * 6.0f / (1.05f * diameter); // m2/m3
float DyeConc               = 0.0033f * 1000.0f / 332.31f;     // mol/m3 

// Physical Constants
const float P1              = 4.59e-10f / DyeConc * BeadConc; //FITC mol / PS area m2
const float P2              = 4.59e-10f / DyeConc * BeadConc; //FITC mol / PS area m2
const float Keq             = 1422.0f;
float Kon1                  = 0 * 100.8f * (W / X) * DyeConc * BeadConc; // 0.10f;
float Kon2                  = 0.0f; //324.8f * (W / X) * DyeConc * BeadConc; // 0.10f;
const float Visc            = 0.0010016f;              //Dynamic viscosity of water at 20C in Pa.s
const float Difusion_dye    = 4.9e-10f; // m2/s From 4.9 × 10−6 cm2 s−1 The diffusion coefficient of fluorescein in water at 21.5°C, as calculated from the Wilke-Chang correlation
const float Difusion_beads  = 1.380649e-23f * (273.15f + 25.0f) / (3.0f * M_PI * Visc * diameter); // m2/s From kB*T/(3*pi*visc*d) kB=1.380649×10−23 J⋅K−1

// Derived values for Crank-Nicolson implicit method
const float restime       = W * H * L / Flow; // seconds
const float dt            = restime / Z; // seconds
const float T_dye         = Difusion_dye * restime / (W * W);
const float T_beads       = Difusion_beads * restime / (W * W);
const float dT_dye        = Difusion_dye * dt / (W * W); // T = Dt/l^2
const float dT_beads      = Difusion_beads * dt / (W * W);
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
    float* DD{new float[M]{}};
    float* ND{new float[M]{}};
    float* NDD{new float[M]{}};
    float* NS{new float[M]{}};

    //initialize the 3 Concentration Arrays.
    for (int z = 0; z < Z; z++) {
        SetInletCondit(D, DD, ND, NS, NDD, z);
    }

    Eigen::MatrixXf m(X,X);
    m.setZero();
    Eigen::MatrixXf m_beads(X,X);
    m_beads.setZero();
    Eigen::MatrixXf E(X,1);
    E.setZero();
    Eigen::MatrixXf E_DD(X,1);
    E_DD.setZero();
    Eigen::MatrixXf E_NS(X,1);
    E_NS.setZero();
    Eigen::MatrixXf E_ND(X,1);
    E_ND.setZero();
    Eigen::MatrixXf E_NDD(X,1);
    E_NDD.setZero();
    Eigen::MatrixXf Solution(X,1);
    Solution.setZero();
    Eigen::MatrixXf Solution_DD(X,1);
    Solution_DD.setZero();
    Eigen::MatrixXf Solution_NS(X,1);
    Solution_NS.setZero();
    Eigen::MatrixXf Solution_ND(X,1);
    Solution_ND.setZero();
    Eigen::MatrixXf Solution_NDD(X,1);
    Solution_NDD.setZero();

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

    float reactionrate1 = 0.0f;
    float reactionrate2 = 0.0f;

    for (int z = 1; z < Z; z++) {
        for (int i = 0; i < X; i++) {
            reactionrate1 = Kon1 * dt * (D[i + (z - 1) * X] * NS[i + (z - 1) * X] - ND[i + (z - 1) * X] / Keq);
            reactionrate2 = 0.0f; //Kon2 * dt * (DD[i + (z - 1) * X] * NDD[i + (z - 1) * X]);

            if (i == 0) {
                E(i,0)      = (1.00f - r_dye)   * D[i + (z - 1) * X]    + r_dye   * D[i + (z - 1) * X + 1]                                                              - reactionrate1;
                E_DD(i,0)   = (1.00f - r_beads) * DD[i + (z - 1) * X]   + r_beads * DD[i + (z - 1) * X + 1]                                                             + reactionrate1      - reactionrate2;
                E_NS(i,0)   = (1.00f - r_beads) * NS[i + (z - 1) * X]   + r_beads * NS[i + (z - 1) * X + 1]                                                             - reactionrate1 / P1;
                E_NDD(i,0)  = (1.00f - r_beads) * NDD[i + (z - 1) * X]  + r_beads * NDD[i + (z - 1) * X + 1]                                                            + reactionrate1 / P1 - reactionrate2 / P2;
                E_ND(i,0)   = (1.00f - r_beads) * ND[i + (z - 1) * X]   + r_beads * ND[i + (z - 1) * X + 1]                                                             + reactionrate1 / P1;
                
            } else if (i == X - 1) {
                E(i,0)      =                                           r_dye   * D[i + (z - 1) * X - 1]                    + (1.00f - r_dye)   * D[i + (z - 1) * X]    - reactionrate1;
                E_DD(i,0)   =                                           r_beads * DD[i + (z - 1) * X - 1]                   + (1.00f - r_beads) * DD[i + (z - 1) * X]   + reactionrate1 - reactionrate2;
                E_NS(i,0)   =                                           r_beads * NS[i + (z - 1) * X - 1]                   + (1.00f - r_beads) * NS[i + (z - 1) * X]   - reactionrate1 / P1;
                E_NDD(i,0)  =                                           r_beads * NDD[i + (z - 1) * X - 1]                  + (1.00f - r_beads) * NDD[i + (z - 1) * X]  + reactionrate1 / P1 - reactionrate2 / P2;
                E_ND(i,0)   =                                           r_beads * ND[i + (z - 1) * X - 1]                   + (1.00f - r_beads) * ND[i + (z - 1) * X]   + reactionrate1 / P1;
            } else {
                E(i,0)      = r_dye   * D[i + (z - 1) * X - 1]          + (2.00f - 2.00f * r_dye)   * D[i + (z - 1) * X]    + r_dye   * D[i + (z - 1) * X + 1]          - reactionrate1;
                E_DD(i,0)   = r_beads * DD[i + (z - 1) * X - 1]         + (2.00f - 2.00f * r_beads) * DD[i + (z - 1) * X]   + r_beads * DD[i + (z - 1) * X + 1]         + reactionrate1 - reactionrate2;
                E_NS(i,0)   = r_beads * NS[i + (z - 1) * X - 1]         + (2.00f - 2.00f * r_beads) * NS[i + (z - 1) * X]   + r_beads * NS[i + (z - 1) * X + 1]         - reactionrate1 / P1;
                E_NDD(i,0)  = r_beads * NDD[i + (z - 1) * X - 1]        + (2.00f - 2.00f * r_beads) * NDD[i + (z - 1) * X]  + r_beads * NDD[i + (z - 1) * X + 1]        + reactionrate1 / P1 - reactionrate2 / P2;
                E_ND(i,0)   = r_beads * ND[i + (z - 1) * X - 1]         + (2.00f - 2.00f * r_beads) * ND[i + (z - 1) * X]   + r_beads * ND[i + (z - 1) * X + 1]         + reactionrate1 / P1;
            }
        }

        Solution        = m*E;
        Solution_DD     = m_beads*E_DD;
        Solution_NS     = m_beads*E_NS;
        Solution_ND     = m_beads*E_ND;
        Solution_NDD    = m_beads*E_NDD;

        for (int i = 0; i < X; i++) {
            D[i + z * X]    = max(min(Solution(i,0), 1.00f), 0.00f);
            NS[i + z * X]   = max(min(Solution_NS(i,0), 1.00f), 0.00f);
            ND[i + z * X]   = max(Solution_ND(i,0), 0.00f);
            NDD[i + z * X]  = max(Solution_NDD(i,0), 0.00f);
            DD[i + z * X]   = ND[i + z * X] * P1;
        }
        
    }
    SaveExcelOutput(Z-1, D, DD, ND, NDD, NS);
    std::cout << "Program finished " << endl;
    return 0;
}

void
SetInletCondit(float* D, float* DD, float* ND, float* NS, float* NDD, int z)
{
    int i = 0;
    for (int x = 0; x < X; x++) {
        i = x + X * z;

        if (x > (X / 2 - 1)) {
            D[i] = 1.0f; // mol/m3
            NS[i] = 0.0f;
        } else {
            D[i] = 0.0f;
            NS[i] = 1.0f;
        }

        DD[i] = 0.0f;
        ND[i] = 0.0f;
        NDD[i] = 0.0f;
    }
}

void
SaveViewable(float* V, int X, int Y, int Z, string Filename)
{
    ofstream fout;
    fout.open(Filename);
    for (int j = 0; j < Y; j++) {
        for (int i = 0; i < X; i++) {
            fout << V[i + j * X + Z * X * Y] << "   ";
        }
        fout << endl;
    }
    fout.close();
}

void
SaveExcelOutput(int z, float* D, float* DD, float* ND, float* NDD, float* NS)
{
    float Sum;
    float* Point;
    ofstream fout;
    fout.open("ExcelProfiles.txt", std::ofstream::out | std::ofstream::app);
    for (int j = 0; j < 5; j++) {
        fout << Kon1 << "  " << Kon2 << "   " << Flow << "  " << DyeConc << "   " <<  BeadConc << "  ";
        if (j == 0) {
            Point = D;
            fout << "D  ";
        }
        if (j == 1) {
            Point = DD;
            fout << "DD  ";
        }
        if (j == 2) {
            Point = NS;
            fout << "NS  ";
        }
        if (j == 3) {
            Point = NDD;
            fout << "NDD  ";
        }
        if (j == 4) {
            Point = ND;
            fout << "ND ";
        }
        for (int x = 0; x < X; x++) {
            Sum = Point[x + z * X];
            fout << Sum << "    ";
        }
        fout << endl;
    }
    fout.close();
}