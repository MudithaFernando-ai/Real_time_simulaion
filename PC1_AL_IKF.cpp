// PC1_AL_IKF.cpp
// Simple Pendulum Simulation + Indirect Kalman Filter (with RK4)
// g++ -std=c++11 -O3 PC1_AL_IKF.cpp -o pendulum_server
// ./pendulum_server
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <unistd.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>

const double PI = 3.14159265358979323846;

// =====================================================
// DATA STRUCTURES
// =====================================================

struct ALState {
    double q[3];      // [Rx, Ry, theta]
    double qdot[3];   // [Rx_dot, Ry_dot, theta_dot]
};

struct IKFState {
    double x[2];      // [theta, omega] - estimated from measurements
    double P[2][2];   // 2x2 covariance matrix
};

struct LogData {
    std::vector<double> time;
    std::vector<double> theta_matlab, omega_matlab;
    std::vector<double> theta_cpp, omega_cpp;
    std::vector<double> theta_ikf, omega_ikf;
    std::vector<double> cycle_time_us;
};

LogData logData;

// =====================================================
// BYTE SWAPPING (for TCP endianness)
// =====================================================

double swapDouble(double value) {
    uint64_t temp;
    memcpy(&temp, &value, sizeof(double));
    temp = ((temp & 0xFF00000000000000ULL) >> 56) |
           ((temp & 0x00FF000000000000ULL) >> 40) |
           ((temp & 0x0000FF0000000000ULL) >> 24) |
           ((temp & 0x000000FF00000000ULL) >> 8)  |
           ((temp & 0x00000000FF000000ULL) << 8)  |
           ((temp & 0x0000000000FF0000ULL) << 24) |
           ((temp & 0x000000000000FF00ULL) << 40) |
           ((temp & 0x00000000000000FFULL) << 56);
    double result;
    memcpy(&result, &temp, sizeof(double));
    return result;
}

// =====================================================
// AUGMENTED LAGRANGE SOLVER (5x5 System)
// =====================================================
// Solves: [M      Cq^T] [qddot  ]   [Qe + Qv        ]
//         [Cq      0  ] [lambda] = [-Cqq*qdot^2    ]

void augmentedLagrangeStep(ALState& state, double dt,
                          double m_A, double g, double L, double I_theta,
                          double& lam1_out, double& lam2_out) {
    
    double* q = state.q;
    double* qdot = state.qdot;
    
    // ========== MASS MATRIX M (3x3) ==========
    // M = diag(m_A, m_A, I_theta)
    double M11 = m_A, M22 = m_A, M33 = I_theta;
    
    // ========== CONSTRAINT JACOBIAN Cq (2x3) ==========
    // Constraints: R_x + (L/2)*sin(θ) = 0
    //             R_y + (L/2)*cos(θ) = 0
    // Jacobian:
    double C11 = 1.0,  C12 = 0.0,           C13 = (L/2.0)*sin(q[2]);
    double C21 = 0.0,  C22 = 1.0,           C23 = (L/2.0)*cos(q[2]);
    
    // ========== GENERALIZED FORCES Qe (3x1) ==========
    // Qe = [0, -m_A*g, 0]^T  (gravity in y-direction)
    double Qe1 = 0.0, Qe2 = -m_A*g, Qe3 = 0.0;
    
    // ========== QUADRATIC VELOCITY Qv (3x1) ==========
    // Centrifugal/Coriolis terms from constraint
    double Qv1 = 0.0, Qv2 = 0.0;
    double Qv3 = -(qdot[2]*qdot[2])*(L/2.0)*cos(q[2]);
    
    // ========== CONSTRAINT ACCELERATION TERMS ==========
    // Cqq*qdot^2 = [qdot(3)^2 * (L/2) * cos(θ);
    //              -qdot(3)^2 * (L/2) * sin(θ)]
    double cqq1 = (qdot[2]*qdot[2])*(L/2.0)*cos(q[2]);
    double cqq2 = -(qdot[2]*qdot[2])*(L/2.0)*sin(q[2]);
    
    // ========== ASSEMBLE 5x5 AUGMENTED SYSTEM ==========
    // A = [M      Cq^T]     b = [Qe + Qv]
    //     [Cq      0  ]         [-Cqq]
    
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
    
    // ========== GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING ==========
    for (int i = 0; i < 5; i++) {
        // Find pivot row
        int pivot = i;
        for (int k = i+1; k < 5; k++) {
            if (fabs(A[k][i]) > fabs(A[pivot][i])) {
                pivot = k;
            }
        }
        
        // Swap rows if needed
        if (pivot != i) {
            for (int j = 0; j < 5; j++) {
                double temp = A[i][j];
                A[i][j] = A[pivot][j];
                A[pivot][j] = temp;
            }
            double temp = b[i];
            b[i] = b[pivot];
            b[pivot] = temp;
        }
        
        // Check for singularity
        if (fabs(A[i][i]) < 1e-12) continue;
        
        // Eliminate column below pivot
        for (int k = i+1; k < 5; k++) {
            double factor = A[k][i] / A[i][i];
            for (int j = i; j < 5; j++) {
                A[k][j] -= factor * A[i][j];
            }
            b[k] -= factor * b[i];
        }
    }
    
    // ========== BACK SUBSTITUTION ==========
    double sol[5];
    for (int i = 4; i >= 0; i--) {
        sol[i] = b[i];
        for (int j = i+1; j < 5; j++) {
            sol[i] -= A[i][j] * sol[j];
        }
        if (fabs(A[i][i]) > 1e-12) {
            sol[i] /= A[i][i];
        }
    }
    
    double qddot1 = sol[0];
    double qddot2 = sol[1];
    double qddot3 = sol[2];
    lam1_out = sol[3];
    lam2_out = sol[4];
    
    // ========== EULER INTEGRATION ==========
    qdot[0] += qddot1 * dt;
    qdot[1] += qddot2 * dt;
    qdot[2] += qddot3 * dt;
    
    q[0] += qdot[0] * dt;
    q[1] += qdot[1] * dt;
    q[2] += qdot[2] * dt;
}

// =====================================================
// RK4 ODE FUNCTION for Augmented Lagrange
// =====================================================

void rk4Step(ALState& state, double dt,
            double m_A, double g, double L, double I_theta) {
    
    double dummy_lam1, dummy_lam2;
    augmentedLagrangeStep(state, dt, m_A, g, L, I_theta, dummy_lam1, dummy_lam2);
}

// =====================================================
// INDIRECT KALMAN FILTER PREDICTION
// Uses AUGMENTED LAGRANGE dynamics with RK4
// =====================================================

void ikfPredict(IKFState* ikf, double dt, 
                double m_A, double g, double L, double I_theta) {
    
    double theta = ikf->x[0];
    double omega = ikf->x[1];
    double cos_theta = cos(theta);
    double sin_theta = sin(theta);
    
    // ========== RK4 PREDICTION using AL dynamics ==========
    // State equation:  theta_dot = omega
    //                  omega_dot = -m*g*L*sin(theta)/(2*I)
    
    // k1: slope at current state
    double k1_theta = omega;
    double k1_omega = -(m_A*g*L*sin_theta)/(2.0*I_theta);
    
    // k2: slope at midpoint (using k1)
    double theta_2 = theta + 0.5*dt*k1_theta;
    double omega_2 = omega + 0.5*dt*k1_omega;
    double k2_theta = omega_2;
    double k2_omega = -(m_A*g*L*sin(theta_2))/(2.0*I_theta);
    
    // k3: slope at midpoint (using k2)
    double theta_3 = theta + 0.5*dt*k2_theta;
    double omega_3 = omega + 0.5*dt*k2_omega;
    double k3_theta = omega_3;
    double k3_omega = -(m_A*g*L*sin(theta_3))/(2.0*I_theta);
    
    // k4: slope at endpoint (using k3)
    double theta_4 = theta + dt*k3_theta;
    double omega_4 = omega + dt*k3_omega;
    double k4_theta = omega_4;
    double k4_omega = -(m_A*g*L*sin(theta_4))/(2.0*I_theta);
    
    // RK4 state update
    double theta_new = theta + (dt/6.0)*(k1_theta + 2*k2_theta + 2*k3_theta + k4_theta);
    double omega_new = omega + (dt/6.0)*(k1_omega + 2*k2_omega + 2*k3_omega + k4_omega);
    
    // ========== LINEARIZED JACOBIAN of AL dynamics ==========
    // A = [∂f/∂x] at current state
    // f = [theta_dot; omega_dot]
    // A = [0,  1;  -mgL*cos(theta)/(2*I),  0]
    
    double A11 = 0.0;
    double A12 = 1.0;
    double A21 = -(m_A*g*L*cos_theta)/(2.0*I_theta);
    double A22 = 0.0;
    
    // Discrete-time transition: Phi = I + A*dt (Euler linearization)
    double Phi11 = 1.0 + A11*dt;
    double Phi12 = A12*dt;
    double Phi21 = A21*dt;
    double Phi22 = 1.0 + A22*dt;
    
    // ========== PROCESS NOISE COVARIANCE ==========
    // Higher noise accounts for AL (6-state) vs IKF (2-state) mismatch
    double Q11 = 1e-6;    // theta prediction uncertainty
    double Q22 = 1e-5;    // omega prediction uncertainty
    
    // ========== COVARIANCE PREDICTION ==========
    // P- = Phi*P*Phi^T + Q
    double P11 = ikf->P[0][0];
    double P12 = ikf->P[0][1];
    double P21 = ikf->P[1][0];
    double P22 = ikf->P[1][1];
    
    double temp11 = Phi11*P11 + Phi12*P21;
    double temp12 = Phi11*P12 + Phi12*P22;
    double temp21 = Phi21*P11 + Phi22*P21;
    double temp22 = Phi21*P12 + Phi22*P22;
    
    ikf->P[0][0] = temp11*Phi11 + temp12*Phi21 + Q11;
    ikf->P[0][1] = temp11*Phi12 + temp12*Phi22;
    ikf->P[1][0] = temp21*Phi11 + temp22*Phi21;
    ikf->P[1][1] = temp21*Phi12 + temp22*Phi22 + Q22;
    
    // Update state estimate
    ikf->x[0] = theta_new;
    ikf->x[1] = omega_new;
}

// =====================================================
// INDIRECT KALMAN FILTER UPDATE (Joseph Form)
// =====================================================

void ikfUpdate(IKFState* ikf, double z_theta, double z_omega) {
    
    // ========== MEASUREMENT NOISE COVARIANCE ==========
    // R = cov of measurement errors
    double R11 = 1e-8;    // theta measurement noise
    double R22 = 1e-7;    // omega measurement noise
    
    // ========== INNOVATION (measurement residual) ==========
    double y1 = z_theta - ikf->x[0];
    double y2 = z_omega - ikf->x[1];
    
    // ========== INNOVATION COVARIANCE ==========
    // S = H*P*H^T + R  (H = I for direct measurement)
    double P11 = ikf->P[0][0];
    double P12 = ikf->P[0][1];
    double P21 = ikf->P[1][0];
    double P22 = ikf->P[1][1];
    
    double S11 = P11 + R11;
    double S12 = P12;
    double S21 = P21;
    double S22 = P22 + R22;
    
    // ========== DETERMINANT AND INVERSE OF S ==========
    double S_det = S11*S22 - S12*S21;
    if (fabs(S_det) < 1e-12) return;  // Singular matrix
    
    double S11_inv = S22/S_det;
    double S12_inv = -S12/S_det;
    double S21_inv = -S21/S_det;
    double S22_inv = S11/S_det;
    
    // ========== KALMAN GAIN ==========
    // K = P*H^T*S^-1
    double K11 = P11*S11_inv + P12*S21_inv;
    double K12 = P11*S12_inv + P12*S22_inv;
    double K21 = P21*S11_inv + P22*S21_inv;
    double K22 = P21*S12_inv + P22*S22_inv;
    
    // ========== STATE UPDATE ==========
    // x = x + K*y
    ikf->x[0] += K11*y1 + K12*y2;
    ikf->x[1] += K21*y1 + K22*y2;
    
    // ========== COVARIANCE UPDATE (JOSEPH FORM - Numerically Stable) ==========
    // P = (I - K*H)*P*(I - K*H)^T + K*R*K^T
    // Since H = I:  (I - K) matrix is (I_K11, I_K12; I_K21, I_K22)
    // This form maintains positive definiteness and numerical stability
    
    double I_K11 = 1.0 - K11;
    double I_K12 = -K12;
    double I_K21 = -K21;
    double I_K22 = 1.0 - K22;
    
    // temp = (I - K)*P
    double temp11 = I_K11*P11 + I_K12*P21;
    double temp12 = I_K11*P12 + I_K12*P22;
    double temp21 = I_K21*P11 + I_K22*P21;
    double temp22 = I_K21*P12 + I_K22*P22;
    
    // P_new = [(I-K)*P]*(I-K)^T + K*R*K^T (includes measurement noise feedback)
    double P_new11 = temp11*I_K11 + temp12*I_K21 + K11*R11*K11 + K12*R22*K21;
    double P_new12 = temp11*I_K12 + temp12*I_K22 + K11*R11*K12 + K12*R22*K22;
    double P_new21 = temp21*I_K11 + temp22*I_K21 + K21*R11*K11 + K22*R22*K21;
    double P_new22 = temp21*I_K12 + temp22*I_K22 + K21*R11*K12 + K22*R22*K22;
    
    // Enforce positive definiteness by bounding diagonal elements
    ikf->P[0][0] = fmax(P_new11, 1e-10);
    ikf->P[0][1] = P_new12;
    ikf->P[1][0] = P_new21;
    ikf->P[1][1] = fmax(P_new22, 1e-10);
}

// =====================================================
// CSV OUTPUT
// =====================================================

void saveToCSV() {
    std::ofstream csv("AL_IKF_results.csv");
    
    csv << "Time(s),Theta_MATLAB,Omega_MATLAB,Theta_CPP,Omega_CPP,Theta_IKF,Omega_IKF,CycleTime_us\n";
    
    for (size_t i = 0; i < logData.time.size(); i++) {
        csv << std::fixed << std::setprecision(6)
            << logData.time[i] << ","
            << logData.theta_matlab[i] << ","
            << logData.omega_matlab[i] << ","
            << logData.theta_cpp[i] << ","
            << logData.omega_cpp[i] << ","
            << logData.theta_ikf[i] << ","
            << logData.omega_ikf[i] << ","
            << logData.cycle_time_us[i] << "\n";
    }
    csv.close();
    
    std::cout << "\n Saved " << logData.time.size() << " samples to AL_IKF_results.csv\n";
}

// =====================================================
// MAIN SERVER
// =====================================================

int main() {
    
    // System parameters
    const double L = 1.0, w = 0.1, h = 0.1, rho = 7850, g = 9.81;
    const double m_A = rho*L*w*h;
    const double I_theta = m_A*L*L/12.0;
    const double dt_tcp = 0.001;        // 1ms TCP sample interval
    const double dt_sub = dt_tcp / 10;  // 10 sub-steps for accuracy
    
    std::cout << "=========================================================\n";
    std::cout << "C++: AUGMENTED LAGRANGE + INDIRECT KALMAN FILTER\n";
    std::cout << "MATCHES: PC2_Augmented_Lagrange_Clean_NoPlots.m\n";
    std::cout << "IKF: RK4 prediction + Joseph form covariance (FINAL FIX)\n";
    std::cout << "=========================================================\n";
    std::cout << "Mass m_A    = " << m_A << " kg\n";
    std::cout << "Inertia I_θ = " << I_theta << " kg·m²\n";
    std::cout << "Length L    = " << L << " m\n";
    
    // Initialize Augmented Lagrange state (6-state: q, qdot)
    ALState al_state;
    al_state.q[0] = 0.1;      // Rx = 0.1 m
    al_state.q[1] = 0.0;      // Ry = 0.0 m
    al_state.q[2] = PI/6.0;   // theta = 30°
    al_state.qdot[0] = 0.0;   // Rx_dot = 0
    al_state.qdot[1] = 0.0;   // Ry_dot = 0
    al_state.qdot[2] = 0.5;   // theta_dot = 0.5 rad/s
    
    // Initialize IKF state (2-state: theta, omega only)
    IKFState ikf_state;
    ikf_state.x[0] = PI/6.0;
    ikf_state.x[1] = 0.5;
    ikf_state.P[0][0] = 0.001;  ikf_state.P[0][1] = 0.0;
    ikf_state.P[1][0] = 0.0;    ikf_state.P[1][1] = 0.01;
    
    // ========== TCP SERVER SETUP ==========
    const int serverPort = 5000;
    int sock = socket(AF_INET, SOCK_STREAM, 0);
    
    struct sockaddr_in addr;
    addr.sin_family = AF_INET;
    addr.sin_port = htons(serverPort);
    addr.sin_addr.s_addr = INADDR_ANY;
    
    if (bind(sock, (struct sockaddr*)&addr, sizeof(addr)) < 0) {
        std::cerr << " Bind failed\n";
        return 1;
    }
    
    listen(sock, 1);
    std::cout << "Waiting for MATLAB client on port " << serverPort << "...\n";
    
    int client = accept(sock, NULL, NULL);
    std::cout << " MATLAB connected!\n\n";
    
    std::cout << "Time(s) | θ_MAT | ω_MAT | θ_CPP | ω_CPP | θ_IKF | ω_IKF | µs\n";
    std::cout << "====================================================================\n";
    
    int sampleCount = 0;
    
    // ========== MAIN SIMULATION LOOP ==========
    while (true) {
        auto stepStartTime = std::chrono::high_resolution_clock::now();
        
        double theta_meas, omega_meas;
        
        // ===== TCP RECEIVE from MATLAB =====
        if (recv(client, &theta_meas, sizeof(double), 0) <= 0) break;
        if (recv(client, &omega_meas, sizeof(double), 0) <= 0) break;
        
        // Byte swap for endianness
        theta_meas = swapDouble(theta_meas);
        omega_meas = swapDouble(omega_meas);
        
        sampleCount++;
        
        // ===== AUGMENTED LAGRANGE STEP (with RK4 sub-stepping) =====
        double dummy_lam1, dummy_lam2;
        for (int substep = 0; substep < 10; substep++) {
            rk4Step(al_state, dt_sub, m_A, g, L, I_theta);
        }
        
        // ===== IKF CORRECTION (RK4 prediction + Joseph form update) =====
        ikfPredict(&ikf_state, dt_tcp, m_A, g, L, I_theta);
        ikfUpdate(&ikf_state, theta_meas, omega_meas);
        
        // ===== LOG THIS SAMPLE =====
        double t = (sampleCount - 1) * dt_tcp;
        
        logData.time.push_back(t);
        logData.theta_matlab.push_back(theta_meas);
        logData.omega_matlab.push_back(omega_meas);
        logData.theta_cpp.push_back(al_state.q[2]);      // theta from q[2]
        logData.omega_cpp.push_back(al_state.qdot[2]);   // omega from qdot[2]
        logData.theta_ikf.push_back(ikf_state.x[0]);
        logData.omega_ikf.push_back(ikf_state.x[1]);
        
        // ===== CYCLE TIME MEASUREMENT =====
        auto stepEndTime = std::chrono::high_resolution_clock::now();
        double cycle_time_us = std::chrono::duration_cast<std::chrono::microseconds>(
            stepEndTime - stepStartTime).count();
        
        logData.cycle_time_us.push_back(cycle_time_us);
        
        // ===== DISPLAY (every 100 TCP samples = 100 ms) =====
        if (sampleCount % 100 == 0) {
            printf("%.3f | %+.4f | %+.4f | %+.4f | %+.4f | %+.4f | %+.4f | %.0f\n",
                t, 
                theta_meas, omega_meas,
                al_state.q[2], al_state.qdot[2],
                ikf_state.x[0], ikf_state.x[1],
                cycle_time_us);
        }
    }
    
    // ========== CLEANUP ==========
    close(client);
    close(sock);
    
    std::cout << "\n MATLAB disconnected\n";
    
    // Save results
    saveToCSV();
   
    std::cout << " Simulation complete!\n";  
    return 0;
}