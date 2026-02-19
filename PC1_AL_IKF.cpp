// PC1_Numerical model
// C++ Augmented Lagrange model + IKF for simple pendulum
// - C++ has its own pendulum model as simulation model
// - MATLAB sends (theta_m, omega_m, alpha_m) at 1 kHz acts as physical system
// - IKF state: x = [theta_ikf, omega_ikf, alpha_ikf]
//for compile
// g++ -std=c++11 -O3 IKF_Num_17_2_2026.cpp -o pendulum_server
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

// Augmented lagrange State : z = [z1, z2, z3] = [Rx, Ry, theta]
struct ALState {
    double z[3];          // z  = [z1, z2, z3]
    double z_dot[3];      // ż = [z1_dot, z2_dot, z3_dot]
    double z_ddot[3];     // z̈ = [z1_ddot, z2_ddot, z3_ddot]
};

// IKFState = error-state filter state x = [z̃_i, ż̃_i, z̈̃_i] [Eq. (17)]
struct IKFState {
    double x[3];          // x_k      (corrected)      [Eq. (25)], x = [dz_i, dz_i_dot, dz_i_ddot]
    double P[3][3];       // P_k      (corrected)      [Eq. (26)]
    double x_pred[3];     // x^-_{k+1} (predicted)     [Eq. (18)]
    double P_pred[3][3];  // P^-_{k+1} (predicted)     [Eq. (19)]
};

struct LogData {
    std::vector<double> time;
    std::vector<double> theta_m, omega_m, alpha_m;
    std::vector<double> theta_cpp, omega_cpp, alpha_cpp;
    std::vector<double> theta_ikf, omega_ikf, alpha_ikf;
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
// Augmented Lagrange
// ======================================================

// M z̈ + C_z^T λ = Q^e + Q^v
// C_z z̈ + C_zz ż^2 = 0
static void solveAL_z_ddot(const double z[3], const double z_dot[3],
                           double m_A, double g, double L, double I_theta,
                           double Q_i_theta,    // Q_i (1DOF) external torque on z3 [Eq. (38)–(40)]
                           double z_ddot_out[3])
{
    // M = diag(m_A, m_A, I_theta)
    double M11 = m_A, M22 = m_A, M33 = I_theta;

    // C_z = [1, 0, (L/2)sin(theta); 0, 1, (L/2)cos(theta)]
    double C11 = 1.0,  C12 = 0.0,             C13 = (L/2.0)*std::sin(z[2]);
    double C21 = 0.0,  C22 = 1.0,             C23 = (L/2.0)*std::cos(z[2]);

    // Q^e = [0; -m_A*g; 0]
    double Qe1 = 0.0;
    double Qe2 = -m_A * g;
    // Qe3 includes the force correction Q_i_theta per Eq. (38)–(40)
    double Qe3 = 0.0 + Q_i_theta; 

    // Q^v_theta = -(ż3^2)*(L/2)*cos(theta)
    double Qv1 = 0.0;
    double Qv2 = 0.0;
    double Qv3 = -(z_dot[2]*z_dot[2])*(L/2.0)*std::cos(z[2]);

    // C_zz ż^2 = [ż3^2*(L/2)*cos(theta); -ż3^2*(L/2)*sin(theta)]
    double Czz1 = (z_dot[2]*z_dot[2])*(L/2.0)*std::cos(z[2]);
    double Czz2 = -(z_dot[2]*z_dot[2])*(L/2.0)*std::sin(z[2]);

    // 5x5 augmented system [ M C_z^T; C_z 0 ] [z̈; λ] = [Q^e+Q^v; -C_zz]
    double A[5][5] = {
        {M11,  0.0,  0.0,  C11, C21},
        {0.0,  M22,  0.0,  C12, C22},
        {0.0,  0.0,  M33,  C13, C23},
        {C11,  C12,  C13,  0.0, 0.0},
        {C21,  C22,  C23,  0.0, 0.0}
    };

    double b[5] = {
        Qe1 + Qv1,
        Qe2 + Qv2,
        Qe3 + Qv3,
        -Czz1,
        -Czz2
    };

    // Gaussian elimination (partial pivoting)
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

    z_ddot_out[0] = sol[0];
    z_ddot_out[1] = sol[1];
    z_ddot_out[2] = sol[2];
}

// ODE: ẏ = [ż; z̈]
static void augmentedLagrangeODE(const double y[6],
                                 double m_A, double g, double L, double I_theta,
                                 double Q_i_theta,
                                 double dydt[6])
{
    double z[3]      = {y[0], y[1], y[2]};
    double z_dot[3]  = {y[3], y[4], y[5]};
    double z_ddot[3];

    solveAL_z_ddot(z, z_dot, m_A, g, L, I_theta, Q_i_theta, z_ddot);

    dydt[0] = z_dot[0];
    dydt[1] = z_dot[1];
    dydt[2] = z_dot[2];
    dydt[3] = z_ddot[0];
    dydt[4] = z_ddot[1];
    dydt[5] = z_ddot[2];
}

// True RK4 step
static void rk4Step(ALState &state, double dt,
                    double m_A, double g, double L, double I_theta,
                    double Q_i_theta)
{
    double y[6] = {
        state.z[0], state.z[1], state.z[2],
        state.z_dot[0], state.z_dot[1], state.z_dot[2]
    };

    double k1[6], k2[6], k3[6], k4[6], yt[6];

    augmentedLagrangeODE(y, m_A, g, L, I_theta, Q_i_theta, k1);

    for (int i = 0; i < 6; ++i) yt[i] = y[i] + 0.5 * dt * k1[i];
    augmentedLagrangeODE(yt, m_A, g, L, I_theta, Q_i_theta, k2);

    for (int i = 0; i < 6; ++i) yt[i] = y[i] + 0.5 * dt * k2[i];
    augmentedLagrangeODE(yt, m_A, g, L, I_theta, Q_i_theta, k3);

    for (int i = 0; i < 6; ++i) yt[i] = y[i] + dt * k3[i];
    augmentedLagrangeODE(yt, m_A, g, L, I_theta, Q_i_theta, k4);

    for (int i = 0; i < 6; ++i) {
        y[i] = y[i] + (dt / 6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    }

    // write back z and ż
    state.z[0] = y[0]; state.z[1] = y[1]; state.z[2] = y[2];
    state.z_dot[0] = y[3]; state.z_dot[1] = y[4]; state.z_dot[2] = y[5];

    // update z̈ for logging
    double z_ddot_now[3];
    double z_now[3]      = {state.z[0], state.z[1], state.z[2]};
    double z_dot_now[3]  = {state.z_dot[0], state.z_dot[1], state.z_dot[2]};
    solveAL_z_ddot(z_now, z_dot_now, m_A, g, L, I_theta, Q_i_theta, z_ddot_now);
    state.z_ddot[0] = z_ddot_now[0];
    state.z_ddot[1] = z_ddot_now[1];
    state.z_ddot[2] = z_ddot_now[2];
}

// ======================================================
// IKF (Indirect  Kalman filter) 
// ======================================================

// Eq. (20): State transition matrix f_x for [dz_i, dz_i_dot, dz_i_ddot]
static void buildFx_fullEq20(double dt,
                             const double z_i_ref[3], const double z_i_dot_ref[3],
                             double m_A, double g, double L, double I_theta,
                             double Q_i_theta,
                             double f_x[3][3])
{
    const double eps = 1e-8;

    double z_ddot0[3], z_ddot_theta[3], z_ddot_omega[3];
    solveAL_z_ddot(z_i_ref, z_i_dot_ref, m_A, g, L, I_theta, Q_i_theta, z_ddot0);

    // Perturb z3 (theta)
    double z_th[3] = {z_i_ref[0], z_i_ref[1], z_i_ref[2] + eps};
    solveAL_z_ddot(z_th, z_i_dot_ref, m_A, g, L, I_theta, Q_i_theta, z_ddot_theta);

    // Perturb ż3 (omega)
    double z_dot_om[3] = {z_i_dot_ref[0], z_i_dot_ref[1], z_i_dot_ref[2] + eps};
    solveAL_z_ddot(z_i_ref, z_dot_om, m_A, g, L, I_theta, Q_i_theta, z_ddot_omega);

    const double z_ddot_alpha0 = z_ddot0[2];
    const double d_z_ddot_d_z      = (z_ddot_theta[2] - z_ddot_alpha0) / eps; // ∂z̈/∂z  [Eq. (20)]
    const double d_z_ddot_d_z_dot  = (z_ddot_omega[2] - z_ddot_alpha0) / eps; // ∂z̈/∂ż [Eq. (20)]

    const double dt2 = dt * dt;

    // Eq. (20) exact structure specialized to 1DOF
    f_x[0][0] = 1.0 + 0.5 * d_z_ddot_d_z     * dt2;
    f_x[0][1] = dt  + 0.5 * d_z_ddot_d_z_dot * dt2;
    f_x[0][2] = 0.5 * dt2;

    // Velocity error row
    f_x[1][0] = d_z_ddot_d_z     * dt;
    f_x[1][1] = 1.0 + d_z_ddot_d_z_dot * dt;
    f_x[1][2] = dt;

    // Acceleration error row
    f_x[2][0] = 0.0;
    f_x[2][1] = 0.0;
    f_x[2][2] = 1.0;
}

// Ξ matrix: only acceleration-level plant noise (1DOF Eq. (21))
static void buildXi(double sigma_z_ddot_i_sq, double Xi[3][3])
{
    // Eq. (21) specialized: Ξ = diag(0, 0, σ^2_{z̈,i})
    Xi[0][0] = 0.0; Xi[0][1] = 0.0; Xi[0][2] = 0.0;
    Xi[1][0] = 0.0; Xi[1][1] = 0.0; Xi[1][2] = 0.0;
    Xi[2][0] = 0.0; Xi[2][1] = 0.0; Xi[2][2] = sigma_z_ddot_i_sq;
}

// Predict step: x^- = 0 [Eq. (18)], P^- = f_x P f_x^T + Ξ [Eq. (19)]
static void ikfPredict(IKFState &ikf, double dt, double sigma_z_ddot_i_sq,
                       const double z_i_ref[3], const double z_i_dot_ref[3],
                       double m_A, double g, double L, double I_theta,
                       double Q_i_theta)
{
    double f_x[3][3], Xi[3][3];
    buildFx_fullEq20(dt, z_i_ref, z_i_dot_ref, m_A, g, L, I_theta, Q_i_theta, f_x); // Eq. (20)
    buildXi(sigma_z_ddot_i_sq, Xi);                                                 // Eq. (21)

    // Eq. (18): x^- = 0
    ikf.x_pred[0] = 0.0;
    ikf.x_pred[1] = 0.0;
    ikf.x_pred[2] = 0.0;

    // Eq. (19): P^- = f_x P f_x^T + Ξ
    double fP[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                fP[i][j] += f_x[i][k] * ikf.P[k][j];

    double f_x_T[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            f_x_T[i][j] = f_x[j][i];

    double Pp[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                Pp[i][j] += fP[i][k] * f_x_T[k][j];

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            ikf.P_pred[i][j] = Pp[i][j] + Xi[i][j];
}

// Update using ONLY z and ż measurements (no z̈ measurement).
// h_x maps which sensors are used [Eq. (27)–(29)].
static void ikfUpdate_z_zdot_only(IKFState &ikf,
                                  double z_i_ref_scalar, double z_i_dot_ref_scalar,
                                  double z_meas_scalar,   double z_dot_meas_scalar,
                                  double sigma_z_meas_sq, double sigma_z_dot_meas_sq)
{
    // y = o - h(z) [Eq. (22)] with o=[z_meas, ż_meas], h=[z_ref, ż_ref]
    double y[2] = { z_meas_scalar - z_i_ref_scalar,
                    z_dot_meas_scalar - z_i_dot_ref_scalar };

    // Measurement Jacobian h_x for [z, ż] only: H = [1 0 0; 0 1 0]
    double H[2][3] = { {1.0, 0.0, 0.0},
                       {0.0, 1.0, 0.0} };

    // Σ (measurement covariance) [Eq. (23)] subset of Eq. (29)
    double Sigma[2][2] = { {sigma_z_meas_sq, 0.0},
                           {0.0,             sigma_z_dot_meas_sq} };

    // S = H P^- H^T + Σ [Eq. (23)]
    double S00 = ikf.P_pred[0][0] + Sigma[0][0];
    double S01 = ikf.P_pred[0][1] + Sigma[0][1];
    double S10 = ikf.P_pred[1][0] + Sigma[1][0];
    double S11 = ikf.P_pred[1][1] + Sigma[1][1];

    double detS = S00*S11 - S01*S10;
    if (std::fabs(detS) < 1e-18) detS = (detS >= 0.0 ? 1e-18 : -1e-18);

    double invS00 =  S11 / detS;
    double invS01 = -S01 / detS;
    double invS10 = -S10 / detS;
    double invS11 =  S00 / detS;

    // K = P^- H^T S^{-1} [Eq. (24)]
    double K[3][2];

    // P^- H^T = [col0 col1] of P_pred
    double PHT0[3] = { ikf.P_pred[0][0], ikf.P_pred[1][0], ikf.P_pred[2][0] };
    double PHT1[3] = { ikf.P_pred[0][1], ikf.P_pred[1][1], ikf.P_pred[2][1] };

    for (int i = 0; i < 3; ++i) {
        K[i][0] = PHT0[i]*invS00 + PHT1[i]*invS10;
        K[i][1] = PHT0[i]*invS01 + PHT1[i]*invS11;
    }

    // x = 0 + K y [Eq. (25)] (x_pred is 0 by Eq. (18))
    ikf.x[0] = ikf.x_pred[0] + K[0][0]*y[0] + K[0][1]*y[1];
    ikf.x[1] = ikf.x_pred[1] + K[1][0]*y[0] + K[1][1]*y[1];
    ikf.x[2] = ikf.x_pred[2] + K[2][0]*y[0] + K[2][1]*y[1];

    // P = (I - K H) P^- [Eq. (26)]
    double I_minus_KH[3][3] = {
        {1.0 - K[0][0],     -K[0][1],        0.0},
        {    -K[1][0], 1.0 - K[1][1],        0.0},
        {    -K[2][0],     -K[2][1],        1.0}
    };

    double P_new[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                P_new[i][j] += I_minus_KH[i][k] * ikf.P_pred[k][j];

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            ikf.P[i][j] = P_new[i][j];
}

// ======================================================
// CSV
// ======================================================

static void saveToCSV() {
    std::ofstream csv("AL_IKF_results_Num.csv");
    csv << "Time,Theta_MATLAB,Omega_MATLAB,Alpha_MATLAB,"
           "Theta_CPP,Omega_CPP,Alpha_CPP,"
           "Theta_IKF,Omega_IKF,Alpha_IKF,CycleTime_us\n";

    for (size_t i = 0; i < logData.time.size(); ++i) {
        csv << std::fixed << std::setprecision(6)
            << logData.time[i] << ","
            << logData.theta_m[i] << ","
            << logData.omega_m[i] << ","
            << logData.alpha_m[i] << ","
            << logData.theta_cpp[i] << ","
            << logData.omega_cpp[i] << ","
            << logData.alpha_cpp[i] << ","
            << logData.theta_ikf[i] << ","
            << logData.omega_ikf[i] << ","
            << logData.alpha_ikf[i] << ","
            << logData.cycle_us[i] << "\n";
    }

    csv.close();
    std::cout << "\nSaved " << logData.time.size()
              << " samples to AL_IKF_results_Num.csv\n";
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



    const double m_A     = rho * L * w * h;
    const double I_theta = m_A * L * L / 12.0;



    const double dt_tcp = 0.001;        // 1 kHz
    const double dt_sub = dt_tcp / 10;  // AL plant substep count = 10

    // 1) PURE MODEL (Uncorrected, for logging/comparison)
    ALState al_pure{};
    al_pure.z[0] = 0.0;  al_pure.z[1] = 0.0;  al_pure.z[2] = PI/3.0;
    al_pure.z_dot[0] = 0.0; al_pure.z_dot[1] = 0.0; al_pure.z_dot[2] = 0.0;
    al_pure.z_ddot[0] = al_pure.z_ddot[1] = al_pure.z_ddot[2] = 0.0;

    // 2) ESTIMATION MODEL (Corrected by IKF, used for feedback)
    ALState al_est{};
    al_est.z[0] = 0.0;  al_est.z[1] = 0.0;  al_est.z[2] = PI/3.0;
    al_est.z_dot[0] = 0.0; al_est.z_dot[1] = 0.0; al_est.z_dot[2] = 0.0;
    al_est.z_ddot[0] = al_est.z_ddot[1] = al_est.z_ddot[2] = 0.0;

    // Init IKF (error state starts at 0)
    IKFState ikf{};
    ikf.x[0] = 0.0; ikf.x[1] = 0.0; ikf.x[2] = 0.0;

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            ikf.P[i][j] = 0.0;

    ikf.P[0][0] = 1e-3;
    ikf.P[1][1] = 1e-3;
    ikf.P[2][2] = 1e-2;

    // FORCE CORRECTION [Eq. (38)–(40)] but in here no force used
    double Q_i_theta_total = 0.0;   // accumulated Q_i on the joint (1DOF)

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

        // MATLAB sends theta_m, omega_m, alpha_m (3 doubles)
        double theta_m, omega_m, alpha_m;
        if (recv(client, &theta_m, sizeof(double), 0) <= 0) break;
        if (recv(client, &omega_m, sizeof(double), 0) <= 0) break;
        if (recv(client, &alpha_m, sizeof(double), 0) <= 0) break;

        theta_m = swapDouble(theta_m);
        omega_m = swapDouble(omega_m);
        alpha_m = swapDouble(alpha_m);

        ++k;
        double t = k * dt_tcp;

        // =================================================
        // 1) Step Pure AL Model 
        // =================================================
        for (int i = 0; i < 10; ++i) {
            rk4Step(al_pure, dt_sub, m_A, g, L, I_theta, 0.0); // No force correction
        }

        // Pure mechanics 
        double theta_pure_cpp = al_pure.z[2];
        double omega_pure_cpp = al_pure.z_dot[2];
        double alpha_pure_cpp = al_pure.z_ddot[2];

        // =================================================
        // 2) Step Estimation Model (Corrected by IKF)
        // =================================================
        for (int i = 0; i < 10; ++i) {
            // Apply accumulated force correction to sustain acceleration prediction
            rk4Step(al_est, dt_sub, m_A, g, L, I_theta, Q_i_theta_total);
        }

        // Reference (mechanics model) states for IKF
        double z_i_ref[3]      = {al_est.z[0],     al_est.z[1],     al_est.z[2]};
        double z_i_dot_ref[3]  = {al_est.z_dot[0], al_est.z_dot[1], al_est.z_dot[2]};
        double z_i_ref_scalar      = al_est.z[2];
        double z_i_dot_ref_scalar  = al_est.z_dot[2];
        double z_i_ddot_ref_scalar = al_est.z_ddot[2];

        // IKF Predict & Update
        const double sigma_z_ddot_i_sq = 0.00001; // σ^2_{z̈,i} plant noise
        ikfPredict(ikf, dt_tcp, sigma_z_ddot_i_sq,
                   z_i_ref, z_i_dot_ref, m_A, g, L, I_theta, Q_i_theta_total);

        const double sigma_z_meas_sq     = 1e-8; // σ^2_z
        const double sigma_z_dot_meas_sq = 1e-7; // σ^2_{ż}
        ikfUpdate_z_zdot_only(ikf,
                              z_i_ref_scalar,     z_i_dot_ref_scalar,
                              theta_m,            omega_m,
                              sigma_z_meas_sq,    sigma_z_dot_meas_sq);

        // ===================================================================
        // 3) Apply Corrections to al_est ONLY (Eqs. (31), (32), (35))
        // ===================================================================

        // x = [dz_i, dz_i_dot, dz_i_ddot]
        double dz_i      = ikf.x[0];
        double dz_i_dot  = ikf.x[1];
        double dz_i_ddot = ikf.x[2];

        double z_i_corr_scalar      = z_i_ref_scalar      + dz_i;      // corrected z_i
        double z_i_dot_corr_scalar  = z_i_dot_ref_scalar  + dz_i_dot;  // corrected ż_i
        double z_i_ddot_corr_scalar = z_i_ddot_ref_scalar + dz_i_ddot; // corrected z̈_i

        // ===================================================================
        // NEW: FORCE CORRECTION [Eq. (38)–(40)]
        // For 1DOF: I_theta * dz̈_i = ΔQ_i_theta
        // ===================================================================
        double Delta_Q_i_theta = I_theta * dz_i_ddot;

        // Persistently update the external force model Q_i
        Q_i_theta_total += Delta_Q_i_theta;

        // Update Estimation Model (independent coordinate)
        al_est.z[2]      = z_i_corr_scalar;
        al_est.z_dot[2]  = z_i_dot_corr_scalar;
        al_est.z_ddot[2] = z_i_ddot_corr_scalar;

        // Reset error states (indirect filter)
        ikf.x[0] = 0.0;
        ikf.x[1] = 0.0;
        ikf.x[2] = 0.0;

        // ===================================================================
        // 4) Logging
        // ===================================================================

        logData.time.push_back(t);

        // MATLAB Input
        logData.theta_m.push_back(theta_m);
        logData.omega_m.push_back(omega_m);
        logData.alpha_m.push_back(alpha_m);

        // CPP (Pure Open-Loop)
        logData.theta_cpp.push_back(theta_pure_cpp);
        logData.omega_cpp.push_back(omega_pure_cpp);
        logData.alpha_cpp.push_back(alpha_pure_cpp);

        // IKF (Corrected Feedback Loop)
        logData.theta_ikf.push_back(z_i_corr_scalar);
        logData.omega_ikf.push_back(z_i_dot_corr_scalar);
        logData.alpha_ikf.push_back(z_i_ddot_corr_scalar);

        auto dt_us = std::chrono::duration_cast<std::chrono::microseconds>(
                         std::chrono::high_resolution_clock::now() - t_start).count();
        logData.cycle_us.push_back(static_cast<double>(dt_us));

        if (k % 1000 == 0) {
            std::printf(
                "t=%.3f: theta_m=%.4f | theta_pure=%.4f | theta_est=%.4f (alpha_est=%.4f) | Q_i_theta_total=%.4f\n",
                t, theta_m, theta_pure_cpp, z_i_corr_scalar, z_i_ddot_corr_scalar, Q_i_theta_total);
        }
    }

    close(client);
    close(sock);

    saveToCSV();
    return 0;
}

