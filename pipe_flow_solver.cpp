#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
using namespace std;

// Control points

    // 1η Bezier: [0, 0.31]
    static const double X0 = 0.0;
    static const double X1 = 0.1;
    static const double X2 = 0.2;
    static const double X3 = 0.21;
    static const double X4 = 0.31;

    static const double Y0 = 0.7;
    static const double Y1 = 0.7;
    static const double Y2 = 0.65;
    static const double Y3 = 0.636;
    static const double Y4 = 0.636;

    // 2η Bezier: [0.31, 1.0]
    static const double X5 = 0.31;
    static const double X6 = 0.41;
    static const double X7 = 0.6;
    static const double X8 = 0.8;
    static const double X9 = 1.0;

    static const double Y5 = 0.636;
    static const double Y6 = 0.636;
    static const double Y7 = 0.65;
    static const double Y8 = 0.7;
    static const double Y9 = 0.7;

    // 1η Bezier καμπύλη (0 -> 0.31)
    inline void nozzle1_x_S_ders(double t, double& x, double& S, double& dxdt, double& dSdt)
    {
        double omt  = 1.0 - t;
        double omt2 = omt * omt;
        double omt3 = omt2 * omt;
        double omt4 = omt2 * omt2;
        double t2   = t * t;
        double t3   = t2 * t;
        double t4   = t2 * t2;

        // x(t)
        x = X0 * omt4
        + 4.0 * X1 * t * omt3
        + 6.0 * X2 * omt2 * t2
        + 4.0 * X3 * omt * t3
        + X4 * t4;

        // S(t) = y(t)
        S = Y0 * omt4
        + 4.0 * Y1 * t * omt3
        + 6.0 * Y2 * t2 * omt2
        + 4.0 * Y3 * t3 * omt
        + Y4 * t4;

        // x'(t)
        dxdt = 4.0 * (
            (X1 - X0) * omt3
        + 3.0 * (X2 - X1) * t * omt2
        + 3.0 * (X3 - X2) * t2 * omt
        + (X4 - X3) * t3
        );

        // S'(t)
        dSdt = 4.0 * (
            (Y1 - Y0) * omt3
        + 3.0 * (Y2 - Y1) * t * omt2
        + 3.0 * (Y3 - Y2) * t2 * omt
        + (Y4 - Y3) * t3
        );
    }

    // 2η Bezier καμπύλη (0.31 -> 1.0)
    inline void nozzle2_x_S_ders(double t, double& x, double& S, double& dxdt, double& dSdt)
    {
        double omt  = 1.0 - t;
        double omt2 = omt * omt;
        double omt3 = omt2 * omt;
        double omt4 = omt2 * omt2;
        double t2   = t * t;
        double t3   = t2 * t;
        double t4   = t2 * t2;

        x = X5 * omt4
        + 4.0 * X6 * t * omt3
        + 6.0 * X7 * omt2 * t2
        + 4.0 * X8 * omt * t3
        + X9 * t4;

        S = Y5 * omt4
        + 4.0 * Y6 * t * omt3
        + 6.0 * Y7 * t2 * omt2
        + 4.0 * Y8 * t3 * omt
        + Y9 * t4;

        dxdt = 4.0 * (
            (X6 - X5) * omt3
        + 3.0 * (X7 - X6) * t * omt2
        + 3.0 * (X8 - X7) * t2 * omt
        + (X9 - X8) * t3
        );

        dSdt = 4.0 * (
            (Y6 - Y5) * omt3
        + 3.0 * (Y7 - Y6) * t * omt2
        + 3.0 * (Y8 - Y7) * t2 * omt
        + (Y9 - Y8) * t3
        );
    }


    
template<typename CurveFunc>
double solve_t_for_x(double x_target, CurveFunc curve)
{
    double tL = 0.0, tR = 1.0;

    for (int it = 0; it < 60; ++it) {   
        double tM = 0.5 * (tL + tR);
        double xM, SM, dxdtM, dSdtM;
        curve(tM, xM, SM, dxdtM, dSdtM);

        if (xM < x_target)
            tL = tM;
        else
            tR = tM;
    }

    return 0.5 * (tL + tR);
}

inline void compute_S_and_dSdx(double x, double& S, double& dSdx)
{
    // Αριστερό ευθύ κομμάτι
    if (x <= 0.0) {
        S    = 0.7*0.7;
        dSdx = 0.0;
        return;
    }

    // 1η Bezier: [0, 0.31]
    if (x > 0.0 && x <= 0.31) {
        double t = solve_t_for_x(x, nozzle1_x_S_ders);

        double xval, Sval, dxdt, dSdt;
        nozzle1_x_S_ders(t, xval, Sval, dxdt, dSdt);

        double r    = Sval;                                   
        double drdx = (fabs(dxdt) > 1e-12) ? (dSdt / dxdt) : 0.0;

        S    = r * r;           
        dSdx = 2.0 * r * drdx; 

        return;
    }

    // 2η Bezier: (0.31, 1.0]
    if (x > 0.31 && x <= 1.0) {
        double t = solve_t_for_x(x, nozzle2_x_S_ders);

        double xval, Sval, dxdt, dSdt;
        nozzle2_x_S_ders(t, xval, Sval, dxdt, dSdt);

        double r    = Sval;                                    
        double drdx = (fabs(dxdt) > 1e-12) ? (dSdt / dxdt) : 0.0;

        S    = r * r;           
        dSdx = 2.0 * r * drdx;  

        return;
    }

    // Δεξί ευθύ κομμάτι: (1.0, 1.5]
    if (x > 1.0 && x <= 1.5) {
        S    = 0.7*0.7;
        dSdx = 0.0;
        return;
    }

    
    if (x < -0.3) {
        S    = 0.7*0.7;
        dSdx = 0.0;
    } else { // x > 1.5
        S    = 0.7*0.7;
        dSdx = 0.0;
    }
}

void precompute_geometry(int Ncells,
                         std::vector<double>& x_center,
                         std::vector<double>& S,
                         std::vector<double>& dSdx)
{
    const double x_min = -0.3;
    const double x_max =  1.5;
    const double L     = x_max - x_min;
    const double dx    = L / Ncells;

    x_center.resize(Ncells);
    S.resize(Ncells);
    dSdx.resize(Ncells);

    for (int i = 0; i < Ncells; ++i) {
        double xc = x_min + (i + 0.5) * dx;
        x_center[i] = xc;

        double Si, dSdxi;
        compute_S_and_dSdx(xc, Si, dSdxi);

        S[i]    = Si;
        dSdx[i] = dSdxi;
    }
}








void compute_RHS(
    double U_curr[][3],         
    double RHS[][3],            
    int Ncells,
    double dx,
    double gamma_gas,
    double Rgas,
    const double S[],           
    const double dSdx[],
    double Tt_in,
    double pt_in,         
    double p_out
     
)


{
    
    static double rho[10000], u[10000], E[10000], p[10000], a[10000]; 
    

    
    for (int i = 0; i < Ncells; ++i) {
        rho[i] = U_curr[i][0];
        u[i]   = U_curr[i][1] / rho[i];
        E[i]   = U_curr[i][2] / rho[i];
        p[i]   = (gamma_gas - 1.0) * rho[i] * (E[i] - 0.5 * u[i] * u[i]);
        a[i]   = std::sqrt(gamma_gas * p[i] / rho[i]);
    }


    
    static double F_face[10001][3];  

    for (int j = 0; j <= Ncells; ++j) {
        double UL[3], UR[3];

        if (j == 0)
{
            
            double rho_int = U_curr[0][0];
            double u_int   = U_curr[0][1] / rho_int;
            double E_int   = U_curr[0][2] / rho_int;
            double p_int   = (gamma_gas - 1.0) * rho_int * (E_int - 0.5*u_int*u_int);
            double a_int   = std::sqrt(gamma_gas * p_int / rho_int);
            double M_int   = u_int / a_int;

           

            // Ισεντροπικές σχέσεις για την είσοδο
            double T_in = Tt_in / (1.0 + 0.5 * (gamma_gas - 1.0) * M_int * M_int);
            double p_in = pt_in * std::pow(T_in / Tt_in, gamma_gas/(gamma_gas - 1.0));
            double rho_in = p_in / (Rgas * T_in);

            double a_in = std::sqrt(gamma_gas * p_in / rho_in);
            double u_in = M_int * a_in;

            double E_in = p_in/((gamma_gas - 1.0)*rho_in) + 0.5*u_in*u_in;

            // Ghost cell 
            UL[0] = rho_in;
            UL[1] = rho_in * u_in;
            UL[2] = rho_in * E_in;

            // Interior side
            UR[0] = U_curr[0][0];
            UR[1] = U_curr[0][1];
            UR[2] = U_curr[0][2];
        }

        else if (j == Ncells) {
            
            UL[0] = U_curr[Ncells-1][0];
            UL[1] = U_curr[Ncells-1][1];
            UL[2] = U_curr[Ncells-1][2];

            
            double rho_R = U_curr[Ncells-1][0];
            double u_R   = U_curr[Ncells-1][1] / rho_R;

           
            double p_R = p_out;     

           
            double E_R = p_R / ((gamma_gas - 1.0) * rho_R) + 0.5 * u_R * u_R;

            
            UR[0] = rho_R;
            UR[1] = rho_R * u_R;
            UR[2] = rho_R * E_R;
        }

        else {
            UL[0] = U_curr[j-1][0]; UL[1] = U_curr[j-1][1]; UL[2] = U_curr[j-1][2];
            UR[0] = U_curr[j][0];   UR[1] = U_curr[j][1];   UR[2] = U_curr[j][2];
        }


        double rhoL, uL, pL, EL, aL;
        double rhoR, uR, pR, ER, aR;
        // -----------------------
            // F+ και F- για UL
            // -----------------------
            rhoL = UL[0];
            uL   = UL[1] / rhoL;
            EL   = UL[2] / rhoL;
            pL   = (gamma_gas - 1.0) * rhoL * (EL - 0.5 * uL * uL);
            aL   = std::sqrt(gamma_gas * pL / rhoL);

            double v2L   = uL * uL;
            double c2L   = aL * aL;
            double H_L   = EL + pL / rhoL;

            double lambda1L = uL;        // u
            double lambda2L = uL + aL;   // u+a
            double lambda3L = uL - aL;   // u-a

            double lambda1Lp = 0.5 * (lambda1L + std::fabs(lambda1L));
            double lambda1Lm = 0.5 * (lambda1L - std::fabs(lambda1L));
            double lambda2Lp = 0.5 * (lambda2L + std::fabs(lambda2L));
            double lambda2Lm = 0.5 * (lambda2L - std::fabs(lambda2L));
            double lambda3Lp = 0.5 * (lambda3L + std::fabs(lambda3L));
            double lambda3Lm = 0.5 * (lambda3L - std::fabs(lambda3L));

            double gm1_over_g = (gamma_gas - 1.0) / gamma_gas;
            double one_over_2g = 1.0 / (2.0 * gamma_gas);
            double c2_over_gm1_L = c2L / (gamma_gas - 1.0);

            double FpL[3] = {0.0, 0.0, 0.0};
            double FmL[3] = {0.0, 0.0, 0.0};

            //λ1+ και λ1-
            {
                double coeff_p = gm1_over_g * rhoL * lambda1Lp;
                double coeff_m = gm1_over_g * rhoL * lambda1Lm;

                FpL[0] += coeff_p * 1.0;
                FpL[1] += coeff_p * uL;
                FpL[2] += coeff_p * (0.5 * v2L);

                FmL[0] += coeff_m * 1.0;
                FmL[1] += coeff_m * uL;
                FmL[2] += coeff_m * (0.5 * v2L);
            }
            //λ2+ και λ2-
            {
                double coeff_p = rhoL * one_over_2g * lambda2Lp;
                double coeff_m = rhoL * one_over_2g * lambda2Lm;

                FpL[0] += coeff_p * 1.0;
                FpL[1] += coeff_p * (uL + aL);
                FpL[2] += coeff_p * (0.5 * v2L + c2_over_gm1_L + aL * uL);

                FmL[0] += coeff_m * 1.0;
                FmL[1] += coeff_m * (uL + aL);
                FmL[2] += coeff_m * (0.5 * v2L + c2_over_gm1_L + aL * uL);
            }
            //λ3+ και λ3-
            {
                double coeff_p = rhoL * one_over_2g * lambda3Lp;
                double coeff_m = rhoL * one_over_2g * lambda3Lm;

                FpL[0] += coeff_p * 1.0;
                FpL[1] += coeff_p * (uL - aL);
                FpL[2] += coeff_p * (0.5 * v2L + c2_over_gm1_L - aL * uL);

                FmL[0] += coeff_m * 1.0;
                FmL[1] += coeff_m * (uL - aL);
                FmL[2] += coeff_m * (0.5 * v2L + c2_over_gm1_L - aL * uL);
            }

            
            // F+ και F- για UR
           
            rhoR = UR[0];
            uR   = UR[1] / rhoR;
            ER   = UR[2] / rhoR;
            pR   = (gamma_gas - 1.0) * rhoR * (ER - 0.5 * uR * uR);
            aR   = std::sqrt(gamma_gas * pR / rhoR);

            double v2R   = uR * uR;
            double c2R   = aR * aR;
            double H_R   = ER + pR / rhoR;

            double lambda1R = uR;
            double lambda2R = uR + aR;
            double lambda3R = uR - aR;

            double lambda1Rp = 0.5 * (lambda1R + std::fabs(lambda1R));
            double lambda1Rm = 0.5 * (lambda1R - std::fabs(lambda1R));
            double lambda2Rp = 0.5 * (lambda2R + std::fabs(lambda2R));
            double lambda2Rm = 0.5 * (lambda2R - std::fabs(lambda2R));
            double lambda3Rp = 0.5 * (lambda3R + std::fabs(lambda3R));
            double lambda3Rm = 0.5 * (lambda3R - std::fabs(lambda3R));

            double c2_over_gm1_R = c2R / (gamma_gas - 1.0);

            double FpR[3] = {0.0, 0.0, 0.0};
            double FmR[3] = {0.0, 0.0, 0.0};

            // λ1+ και λ1-
            {
                double coeff_p = gm1_over_g * rhoR * lambda1Rp;
                double coeff_m = gm1_over_g * rhoR * lambda1Rm;

                FpR[0] += coeff_p * 1.0;
                FpR[1] += coeff_p * uR;
                FpR[2] += coeff_p * (0.5 * v2R);

                FmR[0] += coeff_m * 1.0;
                FmR[1] += coeff_m * uR;
                FmR[2] += coeff_m * (0.5 * v2R);
            }
            // λ2+ και λ2-
            {
                double coeff_p = rhoR * one_over_2g * lambda2Rp;
                double coeff_m = rhoR * one_over_2g * lambda2Rm;

                FpR[0] += coeff_p * 1.0;
                FpR[1] += coeff_p * (uR + aR);
                FpR[2] += coeff_p * (0.5 * v2R + c2_over_gm1_R + aR * uR);

                FmR[0] += coeff_m * 1.0;
                FmR[1] += coeff_m * (uR + aR);
                FmR[2] += coeff_m * (0.5 * v2R + c2_over_gm1_R + aR * uR);
            }
            // λ3+ και λ3-
            {
                double coeff_p = rhoR * one_over_2g * lambda3Rp;
                double coeff_m = rhoR * one_over_2g * lambda3Rm;

                FpR[0] += coeff_p * 1.0;
                FpR[1] += coeff_p * (uR - aR);
                FpR[2] += coeff_p * (0.5 * v2R + c2_over_gm1_R - aR * uR);

                FmR[0] += coeff_m * 1.0;
                FmR[1] += coeff_m * (uR - aR);
                FmR[2] += coeff_m * (0.5 * v2R + c2_over_gm1_R - aR * uR);
            }

            
            // F = F^+(UL) + F^-(UR)
            F_face[j][0] = FpL[0] + FmR[0];
            F_face[j][1] = FpL[1] + FmR[1];
            F_face[j][2] = FpL[2] + FmR[2];
        }
    

    //  RHS 
    for (int i = 0; i < Ncells; ++i) {
        double RHS0 = -(F_face[i+1][0] - F_face[i][0]) / dx;
        double RHS1 = -(F_face[i+1][1] - F_face[i][1]) / dx;
        double RHS2 = -(F_face[i+1][2] - F_face[i][2]) / dx;

        double rho_i = rho[i];
        double u_i   = u[i];
        double E_i   = E[i];
        double p_i   = p[i];

        double Si    = S[i];
        double dSdxi = dSdx[i];

        double geom = -(dSdxi / Si);

        double Q0 = geom * (rho_i * u_i);
        double Q1 = geom * (rho_i * u_i * u_i);
        double Q2 = geom * (u_i * (rho_i * E_i + p_i));

        RHS0 += Q0;
        RHS1 += Q1;
        RHS2 += Q2;

        RHS[i][0] = RHS0;
        RHS[i][1] = RHS1;
        RHS[i][2] = RHS2;
    }
}





int main() {
    const double gamma_gas = 1.4;
    const double Tt_in = 288.0;
    const double pt_in = 1.0332e5;
    const double Rgas  = 287.03;   
    double Mis_out = 0.55; 


    double T_out = Tt_in / (1.0 + 0.5 * (gamma_gas-1) * Mis_out * Mis_out);
    double p_out = pt_in * pow(T_out / Tt_in, gamma_gas / (gamma_gas-1));
    double rho_out = p_out / (Rgas * T_out);
    double a_out   = std::sqrt(gamma_gas * Rgas * T_out);
    double u_out   = Mis_out * a_out;
    double E_out   = p_out / ((gamma_gas-1) * rho_out) + 0.5 * u_out * u_out;

    

  


    //ΠΑΡΑΜΕΤΡΟΙ
    const int    Ncells = 301;      
    const double x_min  = -0.3;
    const double x_max  = 1.5;
    const double dx     = (x_max - x_min) / Ncells;
    const double CFL    = 0.9;
    const int    Nsteps = 150000;     
    


    
    double U[Ncells][3];        
    double U_new[Ncells][3];    
    double F_face[Ncells+1][3]; 

    double rho[Ncells], u[Ncells], p[Ncells], E[Ncells], a[Ncells];
    

    // -----------------------------
    // ΑΡΧΙΚΟΠΟΙΗΣΗ 
    // -----------------------------
    for (int i = 0; i < Ncells; ++i) {
        
        rho[i] = 1.2;
        u[i]   = 70.0;         
        p[i]   = 85325.0;
        a[i]   = std::sqrt(gamma_gas * p[i] / rho[i]);
        E[i]   = p[i] / ((gamma_gas - 1.0) * rho[i]) + 0.5 * u[i] * u[i];

        U[i][0] = rho[i];
        U[i][1] = rho[i] * u[i];
        U[i][2] = rho[i] * E[i];
    }

 
    std::vector<double> x_center, S, dSdx;

    precompute_geometry(Ncells, x_center, S, dSdx);
    

    double K1[Ncells][3], K2[Ncells][3], K3[Ncells][3], K4[Ncells][3];
    double U_stage[Ncells][3];
    double RHS[Ncells][3]; 


    std::ofstream convFile("convergence.txt");


    //ΧΡΟΝΙΚΟ LOOP
    for (int it = 0; it < Nsteps; ++it) {
        

        // CFL
        double maxSpeed = 0.0;

        for (int i = 0; i < Ncells; ++i) {
            double rho_i = U[i][0];
            double u_i   = U[i][1] / rho_i;
            double E_i   = U[i][2] / rho_i;
            double p_i   = (gamma_gas - 1.0) * rho_i * (E_i - 0.5 * u_i * u_i);
            double a_i   = std::sqrt(gamma_gas * p_i / rho_i);  

            double speed_i = std::fabs(u_i) + a_i;               

            if (speed_i > maxSpeed)
                maxSpeed = speed_i;
        }

       
        if (maxSpeed < 1e-8) {
            maxSpeed = 1e-8;
        }

        double dt = CFL * dx / maxSpeed;  

        // K1
        compute_RHS(U,
                    RHS,
                    Ncells,
                    dx,
                    gamma_gas,
                    Rgas,
                    S.data(),
                    dSdx.data(),
                    Tt_in,
                    pt_in,
                    p_out);

        for (int i = 0; i < Ncells; ++i) {
            K1[i][0] = RHS[i][0];
            K1[i][1] = RHS[i][1];
            K1[i][2] = RHS[i][2];
        }

        // K2
        for (int i = 0; i < Ncells; ++i) {
            for (int m = 0; m < 3; ++m) {
                U_stage[i][m] = U[i][m] + 0.5 * dt * K1[i][m];
            }
        }

        compute_RHS(U_stage,
                    RHS,
                    Ncells,
                    dx,
                    gamma_gas,
                    Rgas,
                    S.data(),
                    dSdx.data(),
                    Tt_in,
                    pt_in,
                    p_out);

        for (int i = 0; i < Ncells; ++i) {
            K2[i][0] = RHS[i][0];
            K2[i][1] = RHS[i][1];
            K2[i][2] = RHS[i][2];
        }

        // K3
        for (int i = 0; i < Ncells; ++i) {
            for (int m = 0; m < 3; ++m) {
                U_stage[i][m] = U[i][m] + 0.5 * dt * K2[i][m];
            }
        }

        compute_RHS(U_stage,
                    RHS,
                    Ncells,
                    dx,
                    gamma_gas,
                    Rgas,
                    S.data(),
                    dSdx.data(),
                    Tt_in,
                    pt_in,
                    p_out);

        for (int i = 0; i < Ncells; ++i) {
            K3[i][0] = RHS[i][0];
            K3[i][1] = RHS[i][1];
            K3[i][2] = RHS[i][2];
        }

        // K4
        for (int i = 0; i < Ncells; ++i) {
            for (int m = 0; m < 3; ++m) {
                U_stage[i][m] = U[i][m] + dt * K3[i][m];
            }
        }

        compute_RHS(U_stage,
                    RHS,
                    Ncells,
                    dx,
                    gamma_gas,
                    Rgas,
                    S.data(),
                    dSdx.data(),
                    Tt_in,
                    pt_in,
                    p_out);

        for (int i = 0; i < Ncells; ++i) {
            K4[i][0] = RHS[i][0];
            K4[i][1] = RHS[i][1];
            K4[i][2] = RHS[i][2];
        }

        // U_NEW
        for (int i = 0; i < Ncells; ++i) {
            for (int m = 0; m < 3; ++m) {
                U_new[i][m] = U[i][m]
                            + dt / 6.0 * ( K1[i][m]
                                        + 2.0 * K2[i][m]
                                        + 2.0 * K3[i][m]
                                        + K4[i][m] );
            }
        }
        

        

       
        double maxDiff = 0.0;

        for (int i = 0; i < Ncells; ++i) {
            for (int m = 0; m < 3; ++m) {
                double diff = std::fabs(U_new[i][m] - U[i][m]);
                if (diff > maxDiff)
                    maxDiff = diff;
            }
        }

        // std::cout << "it = " << it << "   maxDiff = " << maxDiff << "\n";


        if (it % 100 == 0) {
        convFile << it << "  " << maxDiff << "\n";
        }

        double maxRes[3] = {0.0, 0.0, 0.0};

        for (int i = 0; i < Ncells; ++i) {
            for (int m = 0; m < 3; ++m) {
                double diff = std::fabs(U_new[i][m] - U[i][m]);
                if (diff > maxRes[m]) maxRes[m] = diff;
            }
        }

        std::cout << "it=" << it
                << "  Res0=" << maxRes[0]
                << "  Res1=" << maxRes[1]
                << "  Res2=" << maxRes[2] << "\n";

     


        
        for (int i = 0; i < Ncells; ++i) {
            U[i][0] = U_new[i][0];
            U[i][1] = U_new[i][1];
            U[i][2] = U_new[i][2];
        }


        //Έλεγχος σύγκλισης
        const double tol = 1e-5;  

        if (maxDiff < tol) {
            std::cout << "Συγκλιση στο it = " << it
                    << " με maxDiff = " << maxDiff << std::endl;
            break;
        }

        

    }


    convFile.close();



     
    std::ofstream fout("geometry.txt");
    for (int i = 0; i < Ncells; ++i) {
        fout << x_center[i] << "  " << S[i] << "\n";
    }
    fout.close();

    
   
    

   


    ofstream fm("mach.txt");

    for (int i = 0; i < Ncells; ++i) {
        double rho_i = U[i][0];
        double u_i   = U[i][1] / rho_i;
        double E_i   = U[i][2] / rho_i;
        double p_i   = (gamma_gas - 1.0) * rho_i * (E_i - 0.5 * u_i * u_i);
        double a_i   = std::sqrt(gamma_gas * p_i / rho_i);
        double M_i   = u_i / a_i;

        double x_i = x_min + (i + 0.5) * dx;

        fm << x_i << " " << M_i << "\n";
    }

    fm.close();

        ofstream frho("rho.txt");
    for (int i = 0; i < Ncells; ++i) {
        double x_i   = x_min + (i + 0.5) * dx;
        double rho_i = U[i][0];
        frho << x_i << " " << rho_i << "\n";
    }

    frho.close();




    return 0;
}
