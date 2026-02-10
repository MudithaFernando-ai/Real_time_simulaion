// PC1 numerical
// C++ Augmented Lagrange model + IKF with NUMERICAL Jacobian
// - C++ has its own Augmented Lagrange model  → theta_cpp, omega_cpp, alpha_cpp
// - MATLAB sends (theta_m, omega_m, alpha_m) at 1 kHz (reference)
// - IKF state: x = [theta_ikf, omega_ikf, alpha_ikf]
//   * Correction only in theta_ikf and omega_ikf using MATLAB theta, omega
//   * alpha_ikf is always computed from theta_ikf and omega_ikf
//   * State‑transition Jacobian F is computed NUMERICALLY (finite differences)

//  compile
//      g++ -std=c++11 -O3 IKF_correct_9_2_2026_N.cpp -o pendulum_server1
//  for run
//      ./pendulum_server1

#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <cstdio>
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

LogData logData;

// ------------------------------------------------------
// Utility
// ------------------------------------------------------

double swapDouble(double value) {
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

// ------------------------------------------------------
// Augmented Lagrange step (same physics as MATLAB code)
// ------------------------------------------------------

void augmentedLagrangeStep(ALState &state, double dt,
                           double m_A, double g, double L, double I_theta,
                           double &lam1_out, double &lam2_out) {
    double *q    = state.q;
    double *qdot = state.qdot;
    double *qddot= state.qddot;

    double M11 = m_A, M22 = m_A, M33 = I_theta;

    double C11 = 1.0,  C12 = 0.0,            C13 = (L/2.0)*std::sin(q[2]);
    double C21 = 0.0,  C22 = 1.0,            C23 = (L/2.0)*std::cos(q[2]);

    // Generalized forces: gravity only, no external torque
    double Qe1 = 0.0;
    double Qe2 = -m_A * g;
    double Qe3 = 0.0;

    // Quadratic velocity terms
    double Qv1 = 0.0;
    double Qv2 = 0.0;
    double Qv3 = -(qdot[2]*qdot[2])*(L/2.0)*std::cos(q[2]);

    // Constraint-second-derivative terms
    double cqq1 = (qdot[2]*qdot[2])*(L/2.0)*std::cos(q[2]);
    double cqq2 = -(qdot[2]*qdot[2])*(L/2.0)*std::sin(q[2]);

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

    // Gaussian elimination with partial pivoting
    for (int i = 0; i < 5; ++i) {
        int pivot = i;
        for (int k = i + 1; k < 5; ++k) {
            if (std::fabs(A[k][i]) > std::fabs(A[pivot][i])) pivot = k;
        }
        if (pivot != i) {
            for (int j = 0; j < 5; ++j) {
                double tmp = A[i][j];
                A[i][j] = A[pivot][j];
                A[pivot][j] = tmp;
            }
            double tmpb = b[i];
            b[i] = b[pivot];
            b[pivot] = tmpb;
        }

        if (std::fabs(A[i][i]) < 1e-12) continue;

        for (int k = i + 1; k < 5; ++k) {
            double factor = A[k][i] / A[i][i];
            for (int j = i; j < 5; ++j)
                A[k][j] -= factor * A[i][j];
            b[k] -= factor * b[i];
        }
    }

    double sol[5];
    for (int i = 4; i >= 0; --i) {
        sol[i] = b[i];
        for (int j = i + 1; j < 5; ++j)
            sol[i] -= A[i][j] * sol[j];
        if (std::fabs(A[i][i]) > 1e-12)
            sol[i] /= A[i][i];
    }

    qddot[0] = sol[0];
    qddot[1] = sol[1];
    qddot[2] = sol[2];
    lam1_out  = sol[3];
    lam2_out  = sol[4];

    qdot[0] += qddot[0] * dt;
    qdot[1] += qddot[1] * dt;
    qdot[2] += qddot[2] * dt;

    q[0] += qdot[0] * dt;
    q[1] += qdot[1] * dt;
    q[2] += qdot[2] * dt;
}

void rk4Step(ALState &state, double dt,
             double m_A, double g, double L, double I_theta) {
    double lam1, lam2;
    augmentedLagrangeStep(state, dt, m_A, g, L, I_theta, lam1, lam2);
}

// ------------------------------------------------------
// Alpha model from IKF theta & omega (Augmented Lagrange)
// ------------------------------------------------------

// For a uniform rod of length L, mass m_A, pivot at one end:
// I_pivot = I_theta + m_A (L/2)^2
// Gravity torque about pivot (sign chosen to match AL result):
//   tau_g = + m_A g (L/2) cos(theta)
// alpha = tau_g / I_pivot
double alphaFromThetaOmega(double theta_ikf, double omega_ikf,
                           double g, double L,
                           double m_A, double I_theta) {
    

    double I_pivot = I_theta + m_A * (L*L) / 4.0;
    double tau_g   =  m_A * g * (L / 2.0) * std::cos(theta_ikf);

    return tau_g / I_pivot;
}

// ------------------------------------------------------
// Nonlinear state transition x_next = f(x)
// x = [theta, omega, alpha]
// alpha is recomputed from theta, omega (no independent dynamics)
// ------------------------------------------------------

void stateTransition(const double x[3],
                     double dt,
                     double g, double L,
                     double m_A, double I_theta,
                     double x_next[3]) {
    double theta = x[0];
    double omega = x[1];

    double alpha_k = alphaFromThetaOmega(theta, omega, g, L, m_A, I_theta);

    double theta_next = theta + omega * dt + 0.5 * alpha_k * dt * dt;
    double omega_next = omega + alpha_k * dt;
    double alpha_next = alphaFromThetaOmega(theta_next, omega_next, g, L, m_A, I_theta);

    x_next[0] = theta_next;
    x_next[1] = omega_next;
    x_next[2] = alpha_next;
}

// ------------------------------------------------------
// Numerical Jacobian F = ∂f/∂x via finite differences
// ------------------------------------------------------

void numericalJacobianF(const double x[3],
                        double dt,
                        double g, double L,
                        double m_A, double I_theta,
                        double F[3][3]) {
    double f_base[3];
    stateTransition(x, dt, g, L, m_A, I_theta, f_base);

    const double eps_base = 1e-6;

    for (int j = 0; j < 3; ++j) {
        double x_pert[3] = {x[0], x[1], x[2]};
        double scale = std::max(1.0, std::fabs(x[j]));
        double eps   = eps_base * scale;

        x_pert[j] += eps;

        double f_pert[3];
        stateTransition(x_pert, dt, g, L, m_A, I_theta, f_pert);

        for (int i = 0; i < 3; ++i) {
            F[i][j] = (f_pert[i] - f_base[i]) / eps;
        }
    }
}

// ------------------------------------------------------
// IKF Predict using numerical F
// ------------------------------------------------------

void ikfPredict(IKFState &ikf, double dt,
                double g, double L,
                double m_A, double I_theta) {
    // Nonlinear propagation
    stateTransition(ikf.x, dt, g, L, m_A, I_theta, ikf.x_pred);

    // Numerical Jacobian F at current x
    double F[3][3];
    numericalJacobianF(ikf.x, dt, g, L, m_A, I_theta, F);

    // Covariance prediction: P = F P F^T + Q
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

// ------------------------------------------------------
// IKF Update: measurement z = [theta_m, omega_m]
// Correction only in theta_ikf & omega_ikf; alpha_ikf recomputed
// ------------------------------------------------------

void ikfUpdate(IKFState &ikf,
               double z_theta, double z_omega,
               double g, double L,
               double m_A, double I_theta) {
    double z[2] = {z_theta, z_omega};

    // H: 2x3, measuring theta and omega
    double H[2][3] = {
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0}
    };

    double R[2][2] = {
        {1e-8, 0.0},
        {0.0,  1e-7}
    };

    // y = z - H x_pred
    double y[2];
    y[0] = z[0] - (H[0][0]*ikf.x_pred[0] + H[0][1]*ikf.x_pred[1] + H[0][2]*ikf.x_pred[2]);
    y[1] = z[1] - (H[1][0]*ikf.x_pred[0] + H[1][1]*ikf.x_pred[1] + H[1][2]*ikf.x_pred[2]);

    // S = H P H^T + R  (2x2)
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
    double invS[2][2];
    invS[0][0] =  S[1][1] / detS;
    invS[0][1] = -S[0][1] / detS;
    invS[1][0] = -S[1][0] / detS;
    invS[1][1] =  S[0][0] / detS;

    // P H^T (3x2)
    double PHT[3][2] = {{0,0},{0,0},{0,0}};
    for (int i = 0; i < 3; ++i) {
        PHT[i][0] = ikf.P[i][0];
        PHT[i][1] = ikf.P[i][1];
    }

    // K = PHT * invS (3x2)
    double K[3][2] = {{0,0},{0,0},{0,0}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 2; ++j)
            for (int k = 0; k < 2; ++k)
                K[i][j] += PHT[i][k] * invS[k][j];

    // x = x_pred + K*y
    double x_new[3];
    for (int i = 0; i < 3; ++i) {
        x_new[i] = ikf.x_pred[i]
                 + K[i][0] * y[0]
                 + K[i][1] * y[1];
    }

    // Enforce alpha_ikf from theta_ikf & omega_ikf (ignore MATLAB alpha)
    double theta_corr = x_new[0];
    double omega_corr = x_new[1];
    double alpha_corr = alphaFromThetaOmega(theta_corr, omega_corr, g, L, m_A, I_theta);
    x_new[2] = alpha_corr;

    for (int i = 0; i < 3; ++i)
        ikf.x[i] = x_new[i];

    // P = (I - K H) P
    double I_KH[3][3] = {
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}
    };

    double KH[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 2; ++k)
                KH[i][j] += K[i][k] * H[k][j];

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            I_KH[i][j] -= KH[i][j];

    double P_new[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                P_new[i][j] += I_KH[i][k] * ikf.P[k][j];

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            ikf.P[i][j] = P_new[i][j];
}

// ------------------------------------------------------
// CSV
// ------------------------------------------------------

void saveToCSV() {
    std::ofstream csv("AL_IKF_results_NUM.csv");
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
              << " samples to AL_IKF_results_NUM.csv\n";
}

// ------------------------------------------------------
// MAIN
// ------------------------------------------------------

int main() {
    // System parameters (match MATLAB)
    const double L   = 1.0;
    const double w   = 0.1;
    const double h   = 0.1;
    const double rho = 7850.0;
    const double g   = 9.81;

    const double m_A     = rho * L * w * h;
    const double I_theta = m_A * L * L / 12.0;

    const double dt_tcp = 0.001;        // 1 kHz from MATLAB
    const double dt_sub = dt_tcp / 10;  // AL substep

    // Initialize AL model state (similar to MATLAB)
    ALState al_state{};
    al_state.q[0]    = 0.1;
    al_state.q[1]    = 0.0;
    al_state.q[2]    = PI / 3.0;
    al_state.qdot[0] = 0.0;
    al_state.qdot[1] = 0.0;
    al_state.qdot[2] = 0.5;

    // Initialize IKF state
    IKFState ikf{};
    ikf.x[0] = PI / 3.0; // theta_ikf
    ikf.x[1] = 0.5;      // omega_ikf
    ikf.x[2] = alphaFromThetaOmega(ikf.x[0], ikf.x[1], g, L, m_A, I_theta);

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            ikf.P[i][j] = 0.0;
    ikf.P[0][0] = 1e-3;
    ikf.P[1][1] = 1e-3;
    ikf.P[2][2] = 1e-2;

    // TCP server
    int serverPort = 5000;
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

    std::cout << "MATLAB connected! AL + IKF (numerical Jacobian)\n\n";

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

        // C++ AL model substeps (independent of MATLAB)
        for (int i = 0; i < 10; ++i) {
            rk4Step(al_state, dt_sub, m_A, g, L, I_theta);
        }

        // IKF: predict + update (correction only in theta_ikf & omega_ikf)
        ikfPredict(ikf, dt_tcp, g, L, m_A, I_theta);
        ikfUpdate(ikf, theta_m, omega_m, g, L, m_A, I_theta);

        // Log everything
        logData.time.push_back(t);

        logData.theta_m.push_back(theta_m);
        logData.omega_m.push_back(omega_m);
        logData.alpha_m.push_back(alpha_m);

        logData.theta_cpp.push_back(al_state.q[2]);
        logData.omega_cpp.push_back(al_state.qdot[2]);
        logData.alpha_cpp.push_back(al_state.qddot[2]);

        logData.theta_ikf.push_back(ikf.x[0]);
        logData.omega_ikf.push_back(ikf.x[1]);
        logData.alpha_ikf.push_back(ikf.x[2]);

        auto dt_us = std::chrono::duration_cast<std::chrono::microseconds>(
                         std::chrono::high_resolution_clock::now() - t_start)
                         .count();
        logData.cycle_us.push_back(static_cast<double>(dt_us));

        if (k % 1000 == 0) {
            std::printf("t=%.3f: "
                        "θm=%.4f, ωm=%.4f, αm=%.4f | "
                        "θcpp=%.4f, ωcpp=%.4f, αcpp=%.4f | "
                        "θikf=%.4f, ωikf=%.4f, αikf=%.4f\n",
                        t,
                        theta_m, omega_m, alpha_m,
                        al_state.q[2], al_state.qdot[2], al_state.qddot[2],
                        ikf.x[0], ikf.x[1], ikf.x[2]);
        }
    }

    close(client);
    close(sock);

    saveToCSV();
    return 0;
}

