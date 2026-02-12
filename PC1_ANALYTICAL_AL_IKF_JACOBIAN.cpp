// PC1_Analatical model
// C++ Augmented Lagrange model + IKF for simple pendulum
// - C++ has its own AL model (similar to MATLAB) → theta_cpp, omega_cpp, alpha_cpp
// - MATLAB sends (theta_m, omega_m, alpha_m) at 1 kHz
// - IKF state: x = [theta_ikf, omega_ikf, alpha_ikf]
//   * Correction only in theta_ikf and omega_ikf using MATLAB theta, omega
//   * alpha_ikf is always computed from theta_ikf and omega_ikf (NOT from MATLAB alpha) from 2/2/2026 feedback
//for complie
// g++ -std=c++11 -O3 IKF_correct_9_2_2026_A.cpp -o pendulum_server
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
#include <unistd.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>

const double PI = 3.14159265358979323846;

// ------------------------------------------------------
// Data structures
// ------------------------------------------------------

struct ALState {
    double q[3];      // [Rx, Ry, theta]
    double qdot[3];   // [Rx_dot, Ry_dot, theta_dot]
    double qddot[3];  // [Rx_ddot, Ry_ddot, theta_ddot]
};

struct IKFState {
    double x[3];      // [theta_ikf, omega_ikf, alpha_ikf]
    double P[3][3];   // covariance
    double x_pred[3]; // predicted state
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
// Augmented Lagrange: solve for qddot given q, qdot
// (matches MATLAB augmentedLagrangeODE internal steps) [file:1]
// ======================================================

static void solveAL_qddot(const double q[3], const double qdot[3],
                          double m_A, double g, double L, double I_theta,
                          double qddot_out[3])
{
    // M = diag(m_A, m_A, I_theta)
    double M11 = m_A, M22 = m_A, M33 = I_theta;

    // Cq = [1, 0, (L/2)sin(theta); 0, 1, (L/2)cos(theta)] [file:1]
    double C11 = 1.0,  C12 = 0.0,            C13 = (L/2.0)*std::sin(q[2]);
    double C21 = 0.0,  C22 = 1.0,            C23 = (L/2.0)*std::cos(q[2]);

    // Qe = [0; -m_A*g; 0] [file:1]
    double Qe1 = 0.0;
    double Qe2 = -m_A * g;
    double Qe3 = 0.0;

    // Qv_theta = -(omega^2)*(L/2)*cos(theta) [file:1]
    double Qv1 = 0.0;
    double Qv2 = 0.0;
    double Qv3 = -(qdot[2]*qdot[2])*(L/2.0)*std::cos(q[2]);

    // Cqq_qdot2 = [omega^2*(L/2)*cos(theta); -omega^2*(L/2)*sin(theta)] [file:1]
    double cqq1 = (qdot[2]*qdot[2])*(L/2.0)*std::cos(q[2]);
    double cqq2 = -(qdot[2]*qdot[2])*(L/2.0)*std::sin(q[2]);

    // 5x5 augmented system [ M Cq^T; Cq 0 ] [qddot; lambda] = [Qe+Qv; -Cqq] [file:1]
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
        -cqq1,
        -cqq2
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

    qddot_out[0] = sol[0];
    qddot_out[1] = sol[1];
    qddot_out[2] = sol[2];
}

// ODE: dydt = [qdot; qddot] (MATLAB augmentedLagrangeODE) [file:1]
static void augmentedLagrangeODE(const double y[6],
                                 double m_A, double g, double L, double I_theta,
                                 double dydt[6])
{
    double q[3]    = {y[0], y[1], y[2]};
    double qdot[3] = {y[3], y[4], y[5]};
    double qddot[3];

    solveAL_qddot(q, qdot, m_A, g, L, I_theta, qddot);

    dydt[0] = qdot[0];
    dydt[1] = qdot[1];
    dydt[2] = qdot[2];
    dydt[3] = qddot[0];
    dydt[4] = qddot[1];
    dydt[5] = qddot[2];
}

// True RK4 step (same structure as MATLAB RK4 loop) [file:1]
static void rk4Step(ALState &state, double dt,
                    double m_A, double g, double L, double I_theta)
{
    double y[6] = {
        state.q[0], state.q[1], state.q[2],
        state.qdot[0], state.qdot[1], state.qdot[2]
    };

    double k1[6], k2[6], k3[6], k4[6], yt[6];

    augmentedLagrangeODE(y, m_A, g, L, I_theta, k1);

    for (int i = 0; i < 6; ++i) yt[i] = y[i] + 0.5 * dt * k1[i];
    augmentedLagrangeODE(yt, m_A, g, L, I_theta, k2);

    for (int i = 0; i < 6; ++i) yt[i] = y[i] + 0.5 * dt * k2[i];
    augmentedLagrangeODE(yt, m_A, g, L, I_theta, k3);

    for (int i = 0; i < 6; ++i) yt[i] = y[i] + dt * k3[i];
    augmentedLagrangeODE(yt, m_A, g, L, I_theta, k4);

    for (int i = 0; i < 6; ++i) {
        y[i] = y[i] + (dt / 6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    }

    // write back q and qdot
    state.q[0] = y[0]; state.q[1] = y[1]; state.q[2] = y[2];
    state.qdot[0] = y[3]; state.qdot[1] = y[4]; state.qdot[2] = y[5];

    // update qddot at the updated state for logging
    double qddot_now[3];
    double q_now[3]    = {state.q[0], state.q[1], state.q[2]};
    double qdot_now[3] = {state.qdot[0], state.qdot[1], state.qdot[2]};
    solveAL_qddot(q_now, qdot_now, m_A, g, L, I_theta, qddot_now);
    state.qddot[0] = qddot_now[0];
    state.qddot[1] = qddot_now[1];
    state.qddot[2] = qddot_now[2];
}

// ======================================================
// IKF: constant-acceleration model + analytical Jacobian
// alpha overwritten after correction by dOmega/dt
// ======================================================

static void ikfPredict(IKFState &ikf, double dt)
{
    const double theta = ikf.x[0];
    const double omega = ikf.x[1];
    const double alpha = ikf.x[2];

    // Prediction (constant acceleration)
    ikf.x_pred[0] = theta + omega * dt + 0.5 * alpha * dt * dt;
    ikf.x_pred[1] = omega + alpha * dt;
    ikf.x_pred[2] = alpha;

    // Analytical Jacobian for this f(x)
    // x_next = [theta + omega*dt + 0.5*alpha*dt^2;
    //           omega + alpha*dt;
    //           alpha]
    double F[3][3] = {
        {1.0, dt, 0.5*dt*dt},
        {0.0, 1.0, dt},
        {0.0, 0.0, 1.0}
    };

    // P = F P F^T + Q
    double P_old[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            P_old[i][j] = ikf.P[i][j];

    double FP[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                FP[i][j] += F[i][k] * P_old[k][j];

    double P_pred[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                P_pred[i][j] += FP[i][k] * F[j][k];

    double Q[3][3] = {
        {1e-8, 0.0,  0.0},
        {0.0,  1e-7, 0.0},
        {0.0,  0.0,  1e-6}
    };

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            ikf.P[i][j] = P_pred[i][j] + Q[i][j];
}

static void ikfUpdate_thetaOmegaOnly(IKFState &ikf, double z_theta, double z_omega)
{
    // Measurement z = [theta_m, omega_m], H = [1 0 0; 0 1 0]
    double H[2][3] = { {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0} };
    double R[2][2] = { {1e-8, 0.0}, {0.0, 1e-7} };

    double y0 = z_theta - (H[0][0]*ikf.x_pred[0] + H[0][1]*ikf.x_pred[1] + H[0][2]*ikf.x_pred[2]);
    double y1 = z_omega - (H[1][0]*ikf.x_pred[0] + H[1][1]*ikf.x_pred[1] + H[1][2]*ikf.x_pred[2]);

    // S = H P H^T + R
    double S[2][2] = {{0,0},{0,0}};
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j) {
            double sum = 0.0;
            for (int k = 0; k < 3; ++k)
                for (int l = 0; l < 3; ++l)
                    sum += H[i][k] * ikf.P[k][l] * H[j][l];
            S[i][j] = sum + R[i][j];
        }

    double detS = S[0][0]*S[1][1] - S[0][1]*S[1][0];
    if (std::fabs(detS) < 1e-18) detS = (detS >= 0 ? 1e-18 : -1e-18);

    double invS[2][2];
    invS[0][0] =  S[1][1] / detS;
    invS[0][1] = -S[0][1] / detS;
    invS[1][0] = -S[1][0] / detS;
    invS[1][1] =  S[0][0] / detS;

    // PHT = P * H^T (H selects first two states)
    double PHT[3][2] = {{0,0},{0,0},{0,0}};
    for (int i = 0; i < 3; ++i) {
        PHT[i][0] = ikf.P[i][0];
        PHT[i][1] = ikf.P[i][1];
    }

    // K = PHT * invS
    double K[3][2] = {{0,0},{0,0},{0,0}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 2; ++j)
            for (int k = 0; k < 2; ++k)
                K[i][j] += PHT[i][k] * invS[k][j];

    // x = x_pred + K*y (but do not correct alpha here)
    double theta_new = ikf.x_pred[0] + K[0][0]*y0 + K[0][1]*y1;
    double omega_new = ikf.x_pred[1] + K[1][0]*y0 + K[1][1]*y1;
    double alpha_new = ikf.x_pred[2];

    ikf.x[0] = theta_new;
    ikf.x[1] = omega_new;
    ikf.x[2] = alpha_new;

    // P = (I - K H) P
    double KH[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 2; ++k)
                KH[i][j] += K[i][k] * H[k][j];

    double I_KH[3][3] = {
        {1.0 - KH[0][0],    -KH[0][1],    -KH[0][2]},
        {   -KH[1][0], 1.0 - KH[1][1],    -KH[1][2]},
        {   -KH[2][0],    -KH[2][1], 1.0 - KH[2][2]}
    };

    double P_new[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                P_new[i][j] += I_KH[i][k] * ikf.P[k][j];

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            ikf.P[i][j] = P_new[i][j];
}

// ======================================================
// CSV
// ======================================================

static void saveToCSV() {
    std::ofstream csv("AL_IKF_results_ANAL_RK4_domega.csv");
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
              << " samples to AL_IKF_results_ANAL_RK4_domega.csv\n";
}

// ======================================================
// MAIN
// ======================================================

int main() {
    // System parameters (match MATLAB) [file:1]
    const double L   = 1.0;
    const double w   = 0.1;
    const double h   = 0.1;
    const double rho = 7850.0;
    const double g   = 9.81;

    const double m_A     = rho * L * w * h;
    const double I_theta = m_A * L * L / 12.0;

    const double dt_tcp = 0.001;        // 1 kHz
    const double dt_sub = dt_tcp / 10;  // AL plant substep count = 10

    // Init C++ AL plant
    ALState al{};
    al.q[0] = 0.1;  al.q[1] = 0.0;  al.q[2] = PI/3.0;
    al.qdot[0] = 0.0; al.qdot[1] = 0.0; al.qdot[2] = 0.5;
    al.qddot[0] = al.qddot[1] = al.qddot[2] = 0.0;

    // Init IKF
    IKFState ikf{};
    ikf.x[0] = PI/3.0;
    ikf.x[1] = 0.5;
    ikf.x[2] = 0.0;   // will be overwritten by dOmega/dt after first update

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            ikf.P[i][j] = 0.0;
    ikf.P[0][0] = 1e-3;
    ikf.P[1][1] = 1e-3;
    ikf.P[2][2] = 1e-2;

    // For alpha_ikf = d(omega_ikf_corrected)/dt
    double omega_prev_corr = 0.0;
    bool have_prev = false;

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

    std::cout << "MATLAB connected! AL(RK4) + IKF(analytical F) + alpha=dOmega/dt\n\n";

    int k = 0;
    while (true) {
        auto t_start = std::chrono::high_resolution_clock::now();

        // MATLAB sends theta_m, omega_m, alpha_m (3 doubles) [file:1]
        double theta_m, omega_m, alpha_m;
        if (recv(client, &theta_m, sizeof(double), 0) <= 0) break;
        if (recv(client, &omega_m, sizeof(double), 0) <= 0) break;
        if (recv(client, &alpha_m, sizeof(double), 0) <= 0) break;

        theta_m = swapDouble(theta_m);
        omega_m = swapDouble(omega_m);
        alpha_m = swapDouble(alpha_m);

        ++k;
        double t = k * dt_tcp;

        // 1) C++ AL plant: RK4 integration (10 substeps)
        for (int i = 0; i < 10; ++i) {
            rk4Step(al, dt_sub, m_A, g, L, I_theta);
        }

        // 2) IKF predict/update (correct theta, omega only)
        ikfPredict(ikf, dt_tcp);
        ikfUpdate_thetaOmegaOnly(ikf, theta_m, omega_m);

        // 3) alpha_ikf from corrected omega derivative
        double omega_corr = ikf.x[1];
        double alpha_ikf = 0.0;
        if (!have_prev) {
            alpha_ikf = 0.0;
            omega_prev_corr = omega_corr;
            have_prev = true;
        } else {
            alpha_ikf = (omega_corr - omega_prev_corr) / dt_tcp;
            omega_prev_corr = omega_corr;
        }
        ikf.x[2] = alpha_ikf;

        // Log
        logData.time.push_back(t);

        logData.theta_m.push_back(theta_m);
        logData.omega_m.push_back(omega_m);
        logData.alpha_m.push_back(alpha_m);

        logData.theta_cpp.push_back(al.q[2]);
        logData.omega_cpp.push_back(al.qdot[2]);
        logData.alpha_cpp.push_back(al.qddot[2]);

        logData.theta_ikf.push_back(ikf.x[0]);
        logData.omega_ikf.push_back(ikf.x[1]);
        logData.alpha_ikf.push_back(ikf.x[2]);

        auto dt_us = std::chrono::duration_cast<std::chrono::microseconds>(
                         std::chrono::high_resolution_clock::now() - t_start).count();
        logData.cycle_us.push_back(static_cast<double>(dt_us));

        if (k % 1000 == 0) {
            std::printf("t=%.3f: θm=%.4f ωm=%.4f αm=%.4f | θikf=%.4f ωikf=%.4f αikf=%.4f | αcpp=%.4f\n",
                        t,
                        theta_m, omega_m, alpha_m,
                        ikf.x[0], ikf.x[1], ikf.x[2],
                        al.qddot[2]);
        }
    }

    close(client);
    close(sock);

    saveToCSV();
    return 0;
}
