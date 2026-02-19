// PC1_Analatical model
// C++ Augmented Lagrange model + IKF for simple pendulum
// - C++ has its own pendulum model as simulation model
// - MATLAB sends (o_z, o_z_dot, o_z_ddot) at 1 kHz acts as physical system
// - IKF state: x = [dz, dz_dot, dz_ddot]
//for compile
// g++ -std=c++11 -O3 Kalman_filter_analytical_17_2_2026.cpp -o pendulum_server
//for run
// ./pendulum_server

#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <algorithm>
#include <cstdint>
#include <unistd.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>

const double PI = 3.14159265358979323846;

// ------------------------------------------------------
// Data structures
// ------------------------------------------------------

struct ALState {
    double Z[3];      // [Rx, Ry, z]
    double Z_dot[3];  // [Rx_dot, Ry_dot, z_dot]
    double Z_ddot[3]; // [Rx_ddot, Ry_ddot, z_ddot]
};

// IKFState = error-state filter state x = [dz, dz_dot, dz_ddot] [Eq. (17)]
struct IKFState {
    double x[3];         // corrected error-state x_k (a posteriori) [Eq. (25)]
    double P[3][3];      // corrected covariance P_k [Eq. (26)]
    double x_minus[3];   // predicted error-state x^-_{k+1} [Eq. (18)]
    double P_minus[3][3];// predicted covariance P^-_{k+1} [Eq. (19)]
};

struct LogData {
    std::vector<double> time;
    std::vector<double> o_z, o_z_dot, o_z_ddot;
    std::vector<double> z_cpp, z_dot_cpp, z_ddot_cpp;
    std::vector<double> z_hat, z_dot_hat, z_ddot_hat;
    std::vector<double> acceleration_error; 
    std::vector<double> cycle_us;
};

static LogData logData;

// ------------------------------------------------------
// Utility: endian swap (MATLAB <-> C++)
// ------------------------------------------------------

static double swapDouble(double value) {
    uint64_t temp;
    std::memcpy(&temp, &value, sizeof(double));
    temp = ((temp & 0xFF00000000000000ULL) >> 56) |
           ((temp & 0x00FF000000000000ULL) >> 40) |
           ((temp & 0x0000FF0000000000ULL) >> 24) |
           ((temp & 0x000000FF00000000ULL) >> 8)  |
           ((temp & 0x00000000FF000000ULL) << 8)  |
           ((temp & 0x0000000000FF0000ULL) << 24) |
           ((temp & 0x000000000000FF00ULL) << 40) |
           ((temp & 0x00000000000000FFULL) << 56);
    double result;
    std::memcpy(&result, &temp, sizeof(double));
    return result;
}

// ======================================================
// Augmented Lagrange: solve for Z_ddot given Z, Z_dot
// (matches MATLAB augmentedLagrangeODE internal steps)
// ======================================================

static void solveAL_Z_ddot(const double Z[3], const double Z_dot[3],
                           double m, double g, double L, double J_zz,
                           double Q_ext, // Added for Eq. 38/40
                           double Z_ddot_out[3])
{
    // M = diag(m, m, J_zz)
    double M11 = m, M22 = m, M33 = J_zz;

    // Phi_Z = [1, 0, (L/2)sin(z); 0, 1, (L/2)cos(z)]
    double Phi_Z11 = 1.0,  Phi_Z12 = 0.0,             Phi_Z13 = (L/2.0)*std::sin(Z[2]);
    double Phi_Z21 = 0.0,  Phi_Z22 = 1.0,             Phi_Z23 = (L/2.0)*std::cos(Z[2]);

    // Qe = [0; -m*g; 0]
    double Qe1 = 0.0;
    double Qe2 = -m * g;
    // Qe3 includes the force correction (Q_ext) per Eq 38/40
    double Qe3 = 0.0 + Q_ext; 

    // Qv_z = -(z_dot^2)*(L/2)*cos(z)
    double Qv1 = 0.0;
    double Qv2 = 0.0;
    double Qv3 = -(Z_dot[2]*Z_dot[2])*(L/2.0)*std::cos(Z[2]);

    // Cqq_Z_dot2 = [z_dot^2*(L/2)*cos(z); -z_dot^2*(L/2)*sin(z)]
    double cqq1 = (Z_dot[2]*Z_dot[2])*(L/2.0)*std::cos(Z[2]);
    double cqq2 = -(Z_dot[2]*Z_dot[2])*(L/2.0)*std::sin(Z[2]);

    // 5x5 augmented system [ M Phi_Z^T; Phi_Z 0 ] [Z_ddot; lambda] = [Qe+Qv; -Cqq]
    double A[5][5] = {
        {M11,  0.0,  0.0,  Phi_Z11, Phi_Z21},
        {0.0,  M22,  0.0,  Phi_Z12, Phi_Z22},
        {0.0,  0.0,  M33,  Phi_Z13, Phi_Z23},
        {Phi_Z11,  Phi_Z12,  Phi_Z13,  0.0, 0.0},
        {Phi_Z21,  Phi_Z22,  Phi_Z23,  0.0, 0.0}
    };

    double b[5] = {
        Qe1 + Qv1,
        Qe2 + Qv2,
        Qe3 + Qv3,
        -cqq1,
        -cqq2
    };

    // Gaussian elimination 
    for (int i = 0; i < 5; ++i) {
        int pivot = i;
        for (int k = i + 1; k < 5; ++k) {
            if (std::fabs(A[k][i]) > std::fabs(A[pivot][i])) pivot = k;
        }
        if (pivot != i) {
            for (int j = 0; j < 5; ++j) std::swap(A[i][j], A[pivot][j]);
            std::swap(b[i], b[pivot]);
        }
        if (std::fabs(A[i][i]) < 1e-12) continue;
        for (int k = i + 1; k < 5; ++k) {
            double factor = A[k][i] / A[i][i];
            for (int j = i; j < 5; ++j) A[k][j] -= factor * A[i][j];
            b[k] -= factor * b[i];
        }
    }

    double sol[5];
    for (int i = 4; i >= 0; --i) {
        sol[i] = b[i];
        for (int j = i + 1; j < 5; ++j) sol[i] -= A[i][j] * sol[j];
        if (std::fabs(A[i][i]) > 1e-12) sol[i] /= A[i][i];
    }

    Z_ddot_out[0] = sol[0];
    Z_ddot_out[1] = sol[1];
    Z_ddot_out[2] = sol[2];
}

// ODE: dydt = [Z_dot; Z_ddot]
static void augmentedLagrangeODE(const double y[6],
                                 double m, double g, double L, double J_zz,
                                 double Q_ext, 
                                 double dydt[6])
{
    double Z[3]    = {y[0], y[1], y[2]};
    double Z_dot[3] = {y[3], y[4], y[5]};
    double Z_ddot[3];

    solveAL_Z_ddot(Z, Z_dot, m, g, L, J_zz, Q_ext, Z_ddot);

    dydt[0] = Z_dot[0];
    dydt[1] = Z_dot[1];
    dydt[2] = Z_dot[2];
    dydt[3] = Z_ddot[0];
    dydt[4] = Z_ddot[1];
    dydt[5] = Z_ddot[2];
}

// True RK4 step
static void rk4Step(ALState &state, double dt,
                    double m, double g, double L, double J_zz,
                    double Q_ext) 
{
    double y[6] = {
        state.Z[0], state.Z[1], state.Z[2],
        state.Z_dot[0], state.Z_dot[1], state.Z_dot[2]
    };

    double k1[6], k2[6], k3[6], k4[6], yt[6];

    augmentedLagrangeODE(y, m, g, L, J_zz, Q_ext, k1);

    for (int i = 0; i < 6; ++i) yt[i] = y[i] + 0.5 * dt * k1[i];
    augmentedLagrangeODE(yt, m, g, L, J_zz, Q_ext, k2);

    for (int i = 0; i < 6; ++i) yt[i] = y[i] + 0.5 * dt * k2[i];
    augmentedLagrangeODE(yt, m, g, L, J_zz, Q_ext, k3);

    for (int i = 0; i < 6; ++i) yt[i] = y[i] + dt * k3[i];
    augmentedLagrangeODE(yt, m, g, L, J_zz, Q_ext, k4);

    for (int i = 0; i < 6; ++i) {
        y[i] = y[i] + (dt / 6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    }

    // write back Z and Z_dot
    state.Z[0] = y[0]; state.Z[1] = y[1]; state.Z[2] = y[2];
    state.Z_dot[0] = y[3]; state.Z_dot[1] = y[4]; state.Z_dot[2] = y[5];

    // update Z_ddot for logging
    double Z_ddot_now[3];
    double Z_now[3]    = {state.Z[0], state.Z[1], state.Z[2]};
    double Z_dot_now[3] = {state.Z_dot[0], state.Z_dot[1], state.Z_dot[2]};
    solveAL_Z_ddot(Z_now, Z_dot_now, m, g, L, J_zz, Q_ext, Z_ddot_now);
    state.Z_ddot[0] = Z_ddot_now[0];
    state.Z_ddot[1] = Z_ddot_now[1];
    state.Z_ddot[2] = Z_ddot_now[2];
}

// ======================================================
// IKF (Indirect / error-state EKF) per paper
// ======================================================

// Eq. (20): State transition matrix for [dz, dz_dot, dz_ddot]

static void buildfx_fullEq20(double dt,
                             const double Z_ref[3], const double Z_dot_ref[3],
                             double m, double g, double L, double J_zz,
                             double Q_ext,
                             double fx[3][3])
{
    // ANALYTICAL JACOBIAN for Single Pendulum
    // Alpha (angular acceleration) equation:
    // J_zz * z_ddot + m*g*(L/2)*sin(z) = Q_ext
    // z_ddot = (Q_ext - m*g*(L/2)*sin(z)) / J_zz

    // dz_ddot_dz = d(z_ddot)/d(z)
    //          = - (m*g*L/2 * cos(z)) / J_zz
    const double dz_ddot_dz = - (m * g * (L/2.0) * std::cos(Z_ref[2])) / J_zz;

    // dz_ddot_dz_dot = d(z_ddot)/d(z_dot)
    //          = 0 (No viscous friction in this model)
    const double dz_ddot_dz_dot = 0.0;

    const double dt2 = dt * dt;

    // Eq. (20) Corrected Structure
    // Row 0: Position error
    fx[0][0] = 1.0 + 0.5 * dz_ddot_dz * dt2;
    fx[0][1] = dt  + 0.5 * dz_ddot_dz_dot   * dt2;
    fx[0][2] = 0.5 * dt2;

    // Row 1: Velocity error
    fx[1][0] = dz_ddot_dz * dt;
    fx[1][1] = 1.0 + dz_ddot_dz_dot * dt;
    fx[1][2] = dt;

    // Row 2: Acceleration error
    fx[2][0] = 0.0;
    fx[2][1] = 0.0;
    fx[2][2] = 1.0;
}

// Sigma_p matrix: only acceleration-level noise (1DOF version of Eq. (21))
static void buildSigma_p(double var_z_ddot, double Sigma_p[3][3])
{
    // Eq. (21) specialized: Sigma_p = diag(0, 0, σ^2_{ẑ,i})
    Sigma_p[0][0] = 0.0; Sigma_p[0][1] = 0.0; Sigma_p[0][2] = 0.0;
    Sigma_p[1][0] = 0.0; Sigma_p[1][1] = 0.0; Sigma_p[1][2] = 0.0;
    Sigma_p[2][0] = 0.0; Sigma_p[2][1] = 0.0; Sigma_p[2][2] = var_z_ddot;
}

// Predict step: x^- = 0 [Eq. (18)], P^- = fx P fx^T + Sigma_p [Eq. (19)]
static void ikfPredict(IKFState &ikf, double dt, double var_z_ddot,
                       const double Z_ref[3], const double Z_dot_ref[3],
                       double m, double g, double L, double J_zz,
                       double Q_ext)
{
    double fx[3][3], Sigma_p[3][3];
    buildfx_fullEq20(dt, Z_ref, Z_dot_ref, m, g, L, J_zz, Q_ext, fx); // Eq. (20)
    buildSigma_p(var_z_ddot, Sigma_p);                                         // Eq. (21)

    // Eq. (18)
    ikf.x_minus[0] = 0.0;
    ikf.x_minus[1] = 0.0;
    ikf.x_minus[2] = 0.0;

    // Eq. (19): P^- = fx * P * fx^T + Sigma_p
    double FP[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                FP[i][j] += fx[i][k] * ikf.P[k][j];

    double fxT[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            fxT[i][j] = fx[j][i];

    double Pp[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                Pp[i][j] += FP[i][k] * fxT[k][j];

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            ikf.P_minus[i][j] = Pp[i][j] + Sigma_p[i][j];
}

// Update using ONLY z and z_dot measurements (no z_ddot measurement).
// hx just maps which sensors you use [Eq. (27)–(29)].
static void ikfUpdate_zZdotOnly(IKFState &ikf,
                                     double z_ref, double z_dot_ref,
                                     double o_z,   double o_z_dot,
                                     double var_z_meas, double var_z_dot_meas)
{
    // Innovation y = o - h(z) [Eq. (22)] with o=[o_z, o_z_dot], h=[z_ref, z_dot_ref]
    double y[2] = { o_z - z_ref, o_z_dot - z_dot_ref };

    // Measurement Jacobian hx for [z, z_dot] only:
    // hx = [1 0 0; 0 1 0]
    double hx[2][3] = { {1.0, 0.0, 0.0},
                        {0.0, 1.0, 0.0} };

    // Sigma_m (measurement covariance) [Eq. (23)] (subset of Eq. (29))
    double Sigma_m[2][2] = { {var_z_meas, 0.0},
                             {0.0, var_z_dot_meas} };

    // S = hx P^- hx^T + Sigma_m [Eq. (23)]
    // Since hx selects first two states, hx P^- hx^T = top-left 2x2 of P_minus
    double S00 = ikf.P_minus[0][0] + Sigma_m[0][0];
    double S01 = ikf.P_minus[0][1] + Sigma_m[0][1];
    double S10 = ikf.P_minus[1][0] + Sigma_m[1][0];
    double S11 = ikf.P_minus[1][1] + Sigma_m[1][1];

    double detS = S00*S11 - S01*S10;
    if (std::fabs(detS) < 1e-18) detS = (detS >= 0.0 ? 1e-18 : -1e-18);

    double invS00 =  S11 / detS;
    double invS01 = -S01 / detS;
    double invS10 = -S10 / detS;
    double invS11 =  S00 / detS;

    // K = P^- hx^T S^{-1} [Eq. (24)]
    // P^- hx^T = [col0 col1] of P_minus
    double K[3][2];

    // First compute PhxT = P_minus * hx^T 
    double PhxT0[3] = { ikf.P_minus[0][0], ikf.P_minus[1][0], ikf.P_minus[2][0] };
    double PhxT1[3] = { ikf.P_minus[0][1], ikf.P_minus[1][1], ikf.P_minus[2][1] };

    // K[:,0] = PhxT[:,0]*invS00 + PhxT[:,1]*invS10
    // K[:,1] = PhxT[:,0]*invS01 + PhxT[:,1]*invS11
    for (int i = 0; i < 3; ++i) {
        K[i][0] = PhxT0[i]*invS00 + PhxT1[i]*invS10;
        K[i][1] = PhxT0[i]*invS01 + PhxT1[i]*invS11;
    }

    // x = 0 + K y [Eq. (25)] (x_minus is 0 by Eq. (18))
    ikf.x[0] = ikf.x_minus[0] + K[0][0]*y[0] + K[0][1]*y[1];
    ikf.x[1] = ikf.x_minus[1] + K[1][0]*y[0] + K[1][1]*y[1];
    ikf.x[2] = ikf.x_minus[2] + K[2][0]*y[0] + K[2][1]*y[1];

    // P = (I - K hx) P^- [Eq. (26)]
    // Compute (I - K hx) explicitly, where hx selects x0,x1
    // Khx = K * hx -> Khx(i,0)=K(i,0), Khx(i,1)=K(i,1), Khx(i,2)=0
    double I_Khx[3][3] = {
        {1.0 - K[0][0],     -K[0][1],        0.0},
        {    -K[1][0], 1.0 - K[1][1],        0.0},
        {    -K[2][0],     -K[2][1],        1.0}
    };

    double P_new[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                P_new[i][j] += I_Khx[i][k] * ikf.P_minus[k][j];

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            ikf.P[i][j] = P_new[i][j];
}

// ======================================================
// CSV
// ======================================================

static void saveToCSV() {
    std::ofstream csv("AL_IKF_results_ANAL.csv");
    // Concatenate strings properly with quotes on each line
    csv << "Time,o_z,o_z_dot,o_z_ddot,"
           "z_CPP,z_dot_CPP,z_ddot_CPP,"
           "z_hat,z_dot_hat,z_ddot_hat,"
           "Acceleration_Error,CycleTime_us\n"; // ADDED COLUMN

    for (size_t i = 0; i < logData.time.size(); ++i) {
        csv << std::fixed << std::setprecision(6)
            << logData.time[i] << ","
            << logData.o_z[i] << ","
            << logData.o_z_dot[i] << ","
            << logData.o_z_ddot[i] << ","
            << logData.z_cpp[i] << ","
            << logData.z_dot_cpp[i] << ","
            << logData.z_ddot_cpp[i] << ","
            << logData.z_hat[i] << ","
            << logData.z_dot_hat[i] << ","
            << logData.z_ddot_hat[i] << ","
            << logData.acceleration_error[i] << "," // WRITE DATA
            << logData.cycle_us[i] << "\n";
    }

    csv.close();
    std::cout << "\nSaved " << logData.time.size()
              << " samples to AL_IKF_results_ANAL.csv\n";
}

// ======================================================
// MAIN
// ======================================================

int main() {
    // System parameters 
    const double L   = 1.0;
    const double w   = 0.1;
    const double h   = 0.1;
    const double rho = 7850.0;
    const double g   = 9.81;

    const double m     = rho * L * w * h;
    const double J_zz = m * L * L / 12.0;

    const double dt_tcp = 0.001;        // 1 kHz
    const double dt_sub = dt_tcp / 10;  // augmented lagrange plant substep count = 10

    // 1) PURE MODEL (Uncorrected, for comparison)
    ALState al_pure{};
    al_pure.Z[0] = 0.0;  al_pure.Z[1] = 0.0;  al_pure.Z[2] = PI/3.0;
    al_pure.Z_dot[0] = 0.0; al_pure.Z_dot[1] = 0.0; al_pure.Z_dot[2] = 0.0;
    al_pure.Z_ddot[0] = al_pure.Z_ddot[1] = al_pure.Z_ddot[2] = 0.0;

    // 2) ESTIMATION MODEL (Corrected by IKF, used for feedback)
    ALState al_est{};
    al_est.Z[0] = 0.0;  al_est.Z[1] = 0.0;  al_est.Z[2] = PI/3.0;
    al_est.Z_dot[0] = 0.0; al_est.Z_dot[1] = 0.0; al_est.Z_dot[2] = 0.0;
    al_est.Z_ddot[0] = al_est.Z_ddot[1] = al_est.Z_ddot[2] = 0.0;

    // Init IKF (error state starts at 0)
    IKFState ikf{};
    ikf.x[0] = 0.0; ikf.x[1] = 0.0; ikf.x[2] = 0.0;

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            ikf.P[i][j] = 0.0;

    ikf.P[0][0] = 1e-3;
    ikf.P[1][1] = 1e-3;
    ikf.P[2][2] = 1e-2;

    // FORCE CORRECTION [Eq 38-40]: Persistent external torque correction
    double Q_ext_total = 0.0;

    // TCP server
    const int serverPort = 5000;
    int sock = socket(AF_INET, SOCK_STREAM, 0);
    if (sock < 0) {
        std::cerr << "Socket creation failed\n";
        return -1;
    }

    sockaddr_in addr{};
    addr.sin_family = AF_INET;
    addr.sin_port   = htons(serverPort);
    addr.sin_addr.s_addr = INADDR_ANY;

    if (bind(sock, (sockaddr*)&addr, sizeof(addr)) < 0) {
        std::cerr << "Bind failed\n";
        close(sock);
        return -1;
    }

    listen(sock, 1);
    std::cout << "Waiting for MATLAB...\n";

    int client = accept(sock, NULL, NULL);
    if (client < 0) {
        std::cerr << "Accept failed\n";
        close(sock);
        return -1;
    }

    std::cout << "MATLAB connected! Dual-State + STRICT Eq(20/35) + FORCE CORRECTION (Eq 38-40)\n\n";

    int k = 0;
    while (true) {
        auto t_start = std::chrono::high_resolution_clock::now();

        // MATLAB sends o_z, o_z_dot, o_z_ddot (3 doubles)
        double o_z, o_z_dot, o_z_ddot;
        if (recv(client, &o_z, sizeof(double), 0) <= 0) break;
        if (recv(client, &o_z_dot, sizeof(double), 0) <= 0) break;
        if (recv(client, &o_z_ddot, sizeof(double), 0) <= 0) break;

        o_z = swapDouble(o_z);
        o_z_dot = swapDouble(o_z_dot);
        o_z_ddot = swapDouble(o_z_ddot);

        ++k;
        double t = k * dt_tcp;

        // =================================================
        // 1) Step Pure AL Model (Unaffected by corrections)
        // =================================================
        for (int i = 0; i < 10; ++i) {
            rk4Step(al_pure, dt_sub, m, g, L, J_zz, 0.0); // No force correction
        }

        // These are the "CPP" variables - purely analytical, no filter
        double z_pure_cpp = al_pure.Z[2];
        double z_dot_pure_cpp = al_pure.Z_dot[2];
        double z_ddot_pure_cpp = al_pure.Z_ddot[2];

        // =================================================
        // 2) Step Estimation Model (Corrected by IKF)
        // =================================================
        for (int i = 0; i < 10; ++i) {
            // Apply accumulated force correction to sustain acceleration prediction
            rk4Step(al_est, dt_sub, m, g, L, J_zz, Q_ext_total);
        }

        // Reference (mechanics model) states for IKF comes from al_est
        double Z_ref[3]    = {al_est.Z[0], al_est.Z[1], al_est.Z[2]};
        double Z_dot_ref[3] = {al_est.Z_dot[0], al_est.Z_dot[1], al_est.Z_dot[2]};
        double z_mech   = al_est.Z[2];
        double z_dot_mech  = al_est.Z_dot[2];
        double z_ddot_mech = al_est.Z_ddot[2];

        // IKF Predict & Update
        const double var_z_ddot = 0.00001; // plant noise
        ikfPredict(ikf, dt_tcp, var_z_ddot, Z_ref, Z_dot_ref, m, g, L, J_zz, Q_ext_total);

        const double var_z_meas = 1e-8;
        const double var_z_dot_meas = 1e-7;
        ikfUpdate_zZdotOnly(ikf, z_mech, z_dot_mech, o_z, o_z_dot,
                                 var_z_meas, var_z_dot_meas);

        // ===================================================================
        // 3) Apply Corrections to al_est ONLY
        // ===================================================================

        // Correct states [Eq. (31), (32)]
        double z_hat = z_mech + ikf.x[0];
        double z_dot_hat = z_dot_mech + ikf.x[1];

        // Correct acceleration [Eq. (35) STRICT]
        // z_ddot_hat = z_ddot_mech + error
        double z_ddot_hat = z_ddot_mech + ikf.x[2];

        // CSV LOGGING: Capture acceleration error BEFORE reset
        double acceleration_error_log = ikf.x[2];

        // ===================================================================
        // NEW: FORCE CORRECTION [Eq 38-40]
        // Calculate torque correction required to produce acceleration error
        // M * delta_z_ddot = delta_Q
        // For 1DOF pendulum: J_zz * dz_ddot = torque_correction
        // ===================================================================
        double dz_ddot = ikf.x[2];
        double delta_Q = J_zz * dz_ddot;

        // Persistently update the external force model
        Q_ext_total += delta_Q;

        // Update Estimation Model
        al_est.Z[2] = z_hat;
        al_est.Z_dot[2] = z_dot_hat;
        al_est.Z_ddot[2] = z_ddot_hat;

        // Reset error states (indirect filter)
        ikf.x[0] = 0.0;
        ikf.x[1] = 0.0;
        ikf.x[2] = 0.0;

        // ===================================================================
        // 4) Logging
        // ===================================================================

        logData.time.push_back(t);

        // MATLAB Input
        logData.o_z.push_back(o_z);
        logData.o_z_dot.push_back(o_z_dot);
        logData.o_z_ddot.push_back(o_z_ddot);

        // CPP (Pure Open-Loop)
        logData.z_cpp.push_back(z_pure_cpp);
        logData.z_dot_cpp.push_back(z_dot_pure_cpp);
        logData.z_ddot_cpp.push_back(z_ddot_pure_cpp);

        // IKF (Corrected Feedback Loop)
        logData.z_hat.push_back(z_hat);
        logData.z_dot_hat.push_back(z_dot_hat);
        logData.z_ddot_hat.push_back(z_ddot_hat);

        // LOG ACCELERATION ERROR
        logData.acceleration_error.push_back(acceleration_error_log);

        auto dt_us = std::chrono::duration_cast<std::chrono::microseconds>(
                         std::chrono::high_resolution_clock::now() - t_start).count();
        logData.cycle_us.push_back(static_cast<double>(dt_us));

        if (k % 1000 == 0) {
            std::printf("t=%.3f: o_z=%.4f | z_pure_cpp=%.4f | z_hat=%.4f (z_ddot_hat=%.4f) | Q_ext=%.4f\n",
                        t, o_z, z_pure_cpp, z_hat, z_ddot_hat, Q_ext_total);
        }
    }

    close(client);
    close(sock);

    saveToCSV();
    return 0;
}
