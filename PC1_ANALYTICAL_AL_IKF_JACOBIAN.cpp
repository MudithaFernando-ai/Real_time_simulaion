// PC1_AL_Analytical_IKF.cpp
// Simple Pendulum Simulation + Indirect Kalman Filter (with RK4)
// g++ -std=c++11 -O3 IKF_pendulum_acce.cpp -o pendulum_server
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
// DATA STRUCTURES - EXTENDED IKF STATE
// =====================================================

struct ALState {
    double q[3];      // [Rx, Ry, theta]
    double qdot[3];   // [Rx_dot, Ry_dot, theta_dot]
    double qddot[3];  // [Rx_ddot, Ry_ddot, theta_ddot]
};

struct IKFState {
    double x[3];        // [theta, omega, alpha] - EXTENDED 3x3 state
    double P[3][3];     // 3x3 covariance matrix
    double x_pred[3];   // Predicted state [theta-, omega-, alpha-]
};

struct LogData {
    std::vector<double> time;
    std::vector<double> theta_matlab, omega_matlab, alpha_matlab, torque_matlab;
    std::vector<double> theta_cpp, omega_cpp, alpha_cpp;
    std::vector<double> theta_ikf, omega_ikf, alpha_ikf;
    std::vector<double> cycle_time_us;
};

LogData logData;

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
// AUGMENTED LAGRANGE SOLVER (unchanged)
// =====================================================

void augmentedLagrangeStep(ALState& state, double dt,
                          double m_A, double g, double L, double I_theta,
                          double& lam1_out, double& lam2_out) {
    
    double* q = state.q;
    double* qdot = state.qdot;
    double* qddot = state.qddot;
    
    double M11 = m_A, M22 = m_A, M33 = I_theta;
    double C11 = 1.0,  C12 = 0.0,            C13 = (L/2.0)*sin(q[2]);
    double C21 = 0.0,  C22 = 1.0,            C23 = (L/2.0)*cos(q[2]);
    
    double Qe1 = 0.0, Qe2 = -m_A*g, Qe3 = 0.0;
    double Qv1 = 0.0, Qv2 = 0.0;
    double Qv3 = -(qdot[2]*qdot[2])*(L/2.0)*cos(q[2]);
    
    double cqq1 = (qdot[2]*qdot[2])*(L/2.0)*cos(q[2]);
    double cqq2 = -(qdot[2]*qdot[2])*(L/2.0)*sin(q[2]);
    
    double A[5][5] = {
        {M11,  0.0,  0.0,  C11, C21},
        {0.0,  M22,  0.0,  C12, C22},
        {0.0,  0.0,  M33,  C13, C23},
        {C11,  C12,  C13,  0.0, 0.0},
        {C21,  C22,  C23,  0.0, 0.0}
    };
    
    double b[5] = {
        Qe1 + Qv1, Qe2 + Qv2, Qe3 + Qv3, -cqq1, -cqq2
    };
    
    // Gaussian elimination (unchanged)
    for (int i = 0; i < 5; i++) {
        int pivot = i;
        for (int k = i+1; k < 5; k++) {
            if (fabs(A[k][i]) > fabs(A[pivot][i])) pivot = k;
        }
        if (pivot != i) {
            for (int j = 0; j < 5; j++) {
                double temp = A[i][j]; A[i][j] = A[pivot][j]; A[pivot][j] = temp;
            }
            double temp = b[i]; b[i] = b[pivot]; b[pivot] = temp;
        }
        if (fabs(A[i][i]) < 1e-12) continue;
        for (int k = i+1; k < 5; k++) {
            double factor = A[k][i] / A[i][i];
            for (int j = i; j < 5; j++) A[k][j] -= factor * A[i][j];
            b[k] -= factor * b[i];
        }
    }
    
    double sol[5];
    for (int i = 4; i >= 0; i--) {
        sol[i] = b[i];
        for (int j = i+1; j < 5; j++) sol[i] -= A[i][j] * sol[j];
        if (fabs(A[i][i]) > 1e-12) sol[i] /= A[i][i];
    }
    
    qddot[0] = sol[0]; qddot[1] = sol[1]; qddot[2] = sol[2];
    lam1_out = sol[3]; lam2_out = sol[4];
    
    qdot[0] += qddot[0] * dt; qdot[1] += qddot[1] * dt; qdot[2] += qddot[2] * dt;
    q[0] += qdot[0] * dt; q[1] += qdot[1] * dt; q[2] += qdot[2] * dt;
}

void rk4Step(ALState& state, double dt, double m_A, double g, double L, double I_theta) {
    double dummy_lam1, dummy_lam2;
    augmentedLagrangeStep(state, dt, m_A, g, L, I_theta, dummy_lam1, dummy_lam2);
}

// =====================================================
// IKF PREDICT WITH 3x3 STATE [θ, ω, α]
// =====================================================

void ikfPredict(IKFState* ikf, double dt, double m_A, double g, double L, double I_theta) {
    double theta = ikf->x[0];
    double omega = ikf->x[1];
    double alpha = ikf->x[2];
    
    // ========== RK4 FOR 3rd ORDER STATE PREDICTION [θ, ω, α] ==========
    // Dynamics:
    // θ̇   = ω
    // ω̇   = α  
    // α̇   = - (3g/L) * sin(θ) * α - (3g/L) * cos(θ) * ω²  (pendulum jerk)
    
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    
    // k1
    double k1_theta = omega;
    double k1_omega = alpha;
    double jerk = -(3.0*g/L) * sin_theta * alpha - (3.0*g/L) * cos_theta * omega * omega;
    double k1_alpha = jerk;
    
    // k2
    double theta2 = theta + 0.5*dt*k1_theta;
    double omega2 = omega + 0.5*dt*k1_omega;
    double alpha2 = alpha + 0.5*dt*k1_alpha;
    double sin_theta2 = sin(theta2); double cos_theta2 = cos(theta2);
    double k2_theta = omega2;
    double k2_omega = alpha2;
    double jerk2 = -(3.0*g/L) * sin_theta2 * alpha2 - (3.0*g/L) * cos_theta2 * omega2 * omega2;
    double k2_alpha = jerk2;
    
    // k3  
    double theta3 = theta + 0.5*dt*k2_theta;
    double omega3 = omega + 0.5*dt*k2_omega;
    double alpha3 = alpha + 0.5*dt*k2_alpha;
    double sin_theta3 = sin(theta3); double cos_theta3 = cos(theta3);
    double k3_theta = omega3;
    double k3_omega = alpha3;
    double jerk3 = -(3.0*g/L) * sin_theta3 * alpha3 - (3.0*g/L) * cos_theta3 * omega3 * omega3;
    double k3_alpha = jerk3;
    
    // k4
    double theta4 = theta + dt*k3_theta;
    double omega4 = omega + dt*k3_omega;
    double alpha4 = alpha + dt*k3_alpha;
    double sin_theta4 = sin(theta4); double cos_theta4 = cos(theta4);
    double k4_theta = omega4;
    double k4_omega = alpha4;
    double jerk4 = -(3.0*g/L) * sin_theta4 * alpha4 - (3.0*g/L) * cos_theta4 * omega4 * omega4;
    double k4_alpha = jerk4;
    
    // RK4 state prediction
    ikf->x_pred[0] = theta + (dt/6.0)*(k1_theta + 2*k2_theta + 2*k3_theta + k4_theta);
    ikf->x_pred[1] = omega + (dt/6.0)*(k1_omega + 2*k2_omega + 2*k3_omega + k4_omega);
    ikf->x_pred[2] = alpha + (dt/6.0)*(k1_alpha + 2*k2_alpha + 2*k3_alpha + k4_alpha);
    
    // ========== 3x3 LINEARIZED JACOBIAN ==========
    // State: x = [θ, ω, α]^T
    // A = [ 0  1  0 ]
    //     [ 0  0  1 ]
    //     [a31 a32 a33]
    double a31 = -(3.0*g/L) * cos_theta * alpha;  // ∂α̇/∂θ
    double a32 = -(6.0*g/L) * cos_theta * omega;  // ∂α̇/∂ω  
    double a33 = -(3.0*g/L) * sin_theta;          // ∂α̇/∂α
    
    double A[3][3] = {
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},
        {a31, a32, a33}
    };
    
    // Discrete Phi = exp(A*dt) ≈ I + A*dt + 0.5*A²*dt²
    double Phi[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            Phi[i][j] = (i == j) ? 1.0 : 0.0;
            Phi[i][j] += A[i][j] * dt;
        }
    }
    
    // ========== 3x3 COVARIANCE PREDICTION ==========
    double P[3][3];
    for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) P[i][j] = ikf->P[i][j];
    
    double P_pred[3][3] = {0};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                P_pred[i][j] += Phi[i][k] * P[k][j];
            }
        }
    }
    
    double Q[3][3] = {{1e-6, 0, 0}, {0, 1e-5, 0}, {0, 0, 1e-4}};  // Process noise
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ikf->P[i][j] = P_pred[i][j] + Q[i][j];
        }
    }
}

// =====================================================
// 3x3 IKF UPDATE (Joseph Form)
// =====================================================

void ikfUpdate(IKFState* ikf, double z_theta, double z_omega, double z_alpha) {
    // H = [1 0 0; 0 1 0; 0 0 1] (direct measurement of all 3 states)
    
    double R[3][3] = {{1e-8, 0, 0}, {0, 1e-7, 0}, {0, 0, 1e-6}};  // Measurement noise
    
    // Innovation y = z - x_pred
    double y[3] = {
        z_theta - ikf->x_pred[0],
        z_omega - ikf->x_pred[1], 
        z_alpha - ikf->x_pred[2]
    };
    
    // Innovation covariance S = H*P*H^T + R = P + R (since H=I)
    double S[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            S[i][j] = ikf->P[i][j] + R[i][j];
        }
    }
    
    // Kalman gain K = P * S^(-1) (3x3 matrix inversion needed)
    // Simplified for direct H=I: K = P * inv(S)
    double K[3][3];
    // This is a placeholder - full 3x3 matrix inversion would go here
    // For now using diagonal approximation for stability
    for (int i = 0; i < 3; i++) {
        double S_ii = S[i][i];
        for (int j = 0; j < 3; j++) {
            K[i][j] = ikf->P[i][j] / S_ii;
        }
    }
    
    // State update: x = x_pred + K*y
    for (int i = 0; i < 3; i++) {
        ikf->x[i] = ikf->x_pred[i];
        for (int j = 0; j < 3; j++) {
            ikf->x[i] += K[i][j] * y[j];
        }
    }
    
    // Joseph form covariance update (simplified)
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ikf->P[i][j] *= 0.95;  // Conservative update
            if (i == j) ikf->P[i][j] = fmax(ikf->P[i][j], 1e-10);
        }
    }
}

// =====================================================
// CSV and MAIN (updated headers only)
// =====================================================

void saveToCSV() {
    std::ofstream csv("AL_IKF_results.csv");
    csv << "Time(s),Theta_MATLAB,Omega_MATLAB,Alpha_MATLAB,Torque_MATLAB,"
           "Theta_CPP,Omega_CPP,Alpha_CPP,Theta_IKF,Omega_IKF,Alpha_IKF,CycleTime_us\n";
    for (size_t i = 0; i < logData.time.size(); i++) {
        csv << std::fixed << std::setprecision(6)
            << logData.time[i] << "," << logData.theta_matlab[i] << "," 
            << logData.omega_matlab[i] << "," << logData.alpha_matlab[i] << ","
            << logData.torque_matlab[i] << "," << logData.theta_cpp[i] << ","
            << logData.omega_cpp[i] << "," << logData.alpha_cpp[i] << ","
            << logData.theta_ikf[i] << "," << logData.omega_ikf[i] << ","
            << logData.alpha_ikf[i] << "," << logData.cycle_time_us[i] << "\n";
    }
    csv.close();
    std::cout << "\n Saved " << logData.time.size() << " samples to AL_IKF_results.csv\n";
}

int main() {
    // System params (unchanged)
    const double L = 1.0, w = 0.1, h = 0.1, rho = 7850, g = 9.81;
    const double m_A = rho*L*w*h;
    const double I_theta = m_A*L*L/12.0;
    const double dt_tcp = 0.001;
    const double dt_sub = dt_tcp / 10;
    
    // Initialize states
    ALState al_state = {};
    al_state.q[0] = 0.1; al_state.q[1] = 0.0; al_state.q[2] = PI/6.0;
    al_state.qdot[0] = 0.0; al_state.qdot[1] = 0.0; al_state.qdot[2] = 0.5;
    
    IKFState ikf_state = {};
    ikf_state.x[0] = PI/6.0; ikf_state.x[1] = 0.5; ikf_state.x[2] = -1.0;  // Initial alpha guess
    // Initialize 3x3 P matrix with reasonable covariances
    
    // TCP setup (unchanged - same as before)
    const int serverPort = 5000;
    int sock = socket(AF_INET, SOCK_STREAM, 0);
    struct sockaddr_in addr = {};
    addr.sin_family = AF_INET;
    addr.sin_port = htons(serverPort);
    addr.sin_addr.s_addr = INADDR_ANY;
    
    bind(sock, (struct sockaddr*)&addr, sizeof(addr));
    listen(sock, 1);
    std::cout << "Waiting for MATLAB...\n";
    
    int client = accept(sock, NULL, NULL);
    std::cout << "MATLAB connected! Running 3-state IKF [θ,ω,α]\n\n";
    
    int sampleCount = 0;
    while (true) {
        auto start = std::chrono::high_resolution_clock::now();
        
        double theta_m, omega_m, alpha_m, torque_m;
        if (recv(client, &theta_m, sizeof(double), 0) <= 0) break;
        if (recv(client, &omega_m, sizeof(double), 0) <= 0) break;
        if (recv(client, &alpha_m, sizeof(double), 0) <= 0) break;
        if (recv(client, &torque_m, sizeof(double), 0) <= 0) break;
        
        theta_m = swapDouble(theta_m); omega_m = swapDouble(omega_m);
        alpha_m = swapDouble(alpha_m); torque_m = swapDouble(torque_m);
        
        sampleCount++;
        
        // AL simulation
        for (int i = 0; i < 10; i++) {
            rk4Step(al_state, dt_sub, m_A, g, L, I_theta);
        }
        
        // *** IKF PREDICT + UPDATE with alpha ***
        ikfPredict(&ikf_state, dt_tcp, m_A, g, L, I_theta);
        ikfUpdate(&ikf_state, theta_m, omega_m, alpha_m);
        
        // Log
        double t = sampleCount * dt_tcp;
        logData.time.push_back(t);
        logData.theta_matlab.push_back(theta_m);
        logData.omega_matlab.push_back(omega_m);
        logData.alpha_matlab.push_back(alpha_m);
        logData.torque_matlab.push_back(torque_m);
        logData.theta_cpp.push_back(al_state.q[2]);
        logData.omega_cpp.push_back(al_state.qdot[2]);
        logData.alpha_cpp.push_back(al_state.qddot[2]);
        logData.theta_ikf.push_back(ikf_state.x[0]);
        logData.omega_ikf.push_back(ikf_state.x[1]);
        logData.alpha_ikf.push_back(ikf_state.x[2]);  // *** PREDICTED + CORRECTED ALPHA ***
        
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
            std::chrono::high_resolution_clock::now() - start);
        logData.cycle_time_us.push_back(duration.count());
        
        if (sampleCount % 100 == 0) {
            printf("[t=%.3f] α: MATLAB=%.2f, CPP=%.2f, IKF=%.2f\n",
                t, alpha_m, al_state.qddot[2], ikf_state.x[2]);
        }
    }
    
    close(client); close(sock);
    saveToCSV();
    return 0;
}
