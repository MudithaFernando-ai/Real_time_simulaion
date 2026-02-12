%% PC2_Augmented_Lagrange_EoM
clear; clc; close all;

%% SYSTEM PARAMETERS
L = 1.0; % Length of body [m]
w = 0.1; % Width [m]
h = 0.1; % Height [m]
rho = 7850; % Density [kg/m^3]
g = 8; % Gravity [m/s^2]

% Mass and inertia
m_A = rho * L * w * h; % Mass m_A
I_theta = m_A * L^2 / 12; % Moment of inertia about COM
fprintf('===============================================\\n');
fprintf('AUGMENTED LAGRANGE - SIMPLE PENDULUM (EoM ALPHA)\\n');
fprintf('===============================================\\n');
fprintf('Length l = %.3f m\\n', L);
fprintf('Mass m_A = %.3f kg\\n', m_A);
fprintf('Inertia I_θ = %.6f kg·m²\\n', I_theta);
fprintf('===============================================\\n\\n');

%% SIMULATION PARAMETERS
dt = 0.001; % Timestep [s] (1 kHz)
t_end = 10.0; % Duration [s]
t = (0:dt:t_end)'; % Time vector
N = length(t);

%% INITIAL CONDITIONS
q0 = [0.1; 0.0; pi/3]; % Initial: small x offset, theta=30°
qdot0 = [0; 0; 0.5]; % Initial velocities
y0 = [q0; qdot0]; % Full state [6x1]

fprintf('Initial Conditions:\\n');
fprintf(' q(0) = [%.3f, %.3f, %.3f] rad\\n', q0);
fprintf(' qdot(0) = [%.3f, %.3f, %.3f] rad/s\\n', qdot0);
fprintf('\\n');

%% TCP/IP CONFIGURATION
serverIP = '169.254.131.136';
serverPort = 5000;

%% CONNECT TO C++ SERVER
try
    tcpipClient = tcpip(serverIP, serverPort, 'NetworkRole', 'client');
    tcpipClient.InputBufferSize = 8192;
    tcpipClient.OutputBufferSize = 8192;
    tcpipClient.Timeout = 30;
    fopen(tcpipClient);
    fprintf(' Connected to C++ Augmented Lagrange server!\\n\\n');
catch ME
    fprintf(' Connection failed: %s\\n', ME.message);
    return;
end

%% SIMULATION STORAGE
y = zeros(6, N); % [R_x, R_y, θ, R_x_dot, R_y_dot, θ_dot]
y(:,1) = y0;
lambda_log = zeros(2, N); % Lagrange multipliers [Fx, Fy]
alpha_log = zeros(1, N); % Angular acceleration (now from EoM)
torque_log = zeros(1, N); % Constraint reaction torque

%% MAIN AUGMENTED LAGRANGE SIMULATION LOOP (MODIFIED ALPHA)
tic;
fprintf('===============================================\\n');
fprintf('AUGMENTED LAGRANGE SIMULATION STARTED (EoM ALPHA)\\n');
fprintf('===============================================\\n');
for k = 1:N-1

    % Current state
    q = y(1:3,k);
    qdot = y(4:6,k);
    
    % ===============================================
    % STEP 1: SIMPLE PENDULUM EQUATION OF MOTION (NEW)
    % ===============================================
    % α = -(m*g*(L/2)*sinθ) / I_θ
    alpha_EoM = -(m_A * g * (L/2) * sin(q(3))) / I_theta;
    
    % ===============================================
    % STEP 2: ASSEMBLE MASS MATRIX M (3x3) - KEPT FOR CONSISTENCY
    % ===============================================
    M = diag([m_A, m_A, I_theta]);
    
    % ===============================================
    % STEP 3: CONSTRAINT JACOBIAN C_q (2x3)
    % ===============================================
    Cq = [1, 0, (L/2)*sin(q(3));
          0, 1, (L/2)*cos(q(3))];
    
    % ===============================================
    % STEP 4: GENERALIZED FORCES Q_e (3x1)
    % ===============================================
    Qe = [0; -m_A*g; 0];
    
    % ===============================================
    % STEP 5: QUADRATIC VELOCITY Q_v (3x1)
    % ===============================================
    Qv_theta = -(qdot(3)^2)*(L/2)*cos(q(3));
    Qv = [0; 0; Qv_theta];
    
    % ===============================================
    % STEP 6: CONSTRAINT ACCELERATION TERMS
    % ===============================================
    Cqq_qdot2 = [qdot(3)^2*(L/2)*cos(q(3));
                 -qdot(3)^2*(L/2)*sin(q(3))];
    
    % ===============================================
    % STEP 7: SOLVE 5x5 AUGMENTED SYSTEM (for lambda, torque)
    % ===============================================
    A = [M, Cq';
         Cq, zeros(2)];
    b = [Qe + Qv; -Cqq_qdot2];
    sol = A \ b;
    
    qddot_full = sol(1:3);  % Full accelerations (for reference)
    lambda = sol(4:5);      % Constraint forces [Fx, Fy] at pivot
    
    % OVERRIDE: Use EoM alpha instead of augmented solve
    qddot = [qddot_full(1); qddot_full(2); alpha_EoM];
    
    % ===============================================
    % COMPUTE REACTION TORQUE about COM
    % ===============================================
    rx = -(L/2)*sin(q(3)); % x-component from pivot to COM
    ry = -(L/2)*cos(q(3)); % y-component from pivot to COM
    torque = rx*lambda(2) - ry*lambda(1); % τ = rx*Fy - ry*Fx
    
    % Store values
    alpha = alpha_EoM;  % NOW FROM EQUATION OF MOTION!
    lambda_log(:,k) = lambda;
    alpha_log(k) = alpha;
    torque_log(k) = torque;
    
    % ===============================================
    % STEP 8: RK4 INTEGRATION (using modified qddot)
    % ===============================================
    k1 = augmentedLagrangeODE_EoM(y(:,k), m_A, g, L, I_theta);
    k2 = augmentedLagrangeODE_EoM(y(:,k) + dt/2*k1, m_A, g, L, I_theta);
    k3 = augmentedLagrangeODE_EoM(y(:,k) + dt/2*k2, m_A, g, L, I_theta);
    k4 = augmentedLagrangeODE_EoM(y(:,k) + dt*k3, m_A, g, L, I_theta);
    
    y_next = y(:,k) + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    y(:,k+1) = y_next;
    
    % ===============================================
    % STEP 9: TCP TRANSMISSION (θ, ω, α_EoM, τ)
    % ===============================================
    theta_current = y_next(3);
    omega_current = y_next(6);
    alpha_current = alpha_EoM;
    torque_current = torque;
    
    try
        fwrite(tcpipClient, theta_current, 'double'); % Angle [rad]
        fwrite(tcpipClient, omega_current, 'double'); % Angular vel [rad/s]
        fwrite(tcpipClient, alpha_current, 'double'); % Angular acc FROM EoM [rad/s²]
        fwrite(tcpipClient, torque_current, 'double'); % Torque [N·m]
        
        if mod(k, 1000) == 0 % Progress every 1s
            fprintf('Step %d: θ=%.3f, ω=%.3f, α_EoM=%.3f, τ=%.3f\\n', ...
                k, theta_current, omega_current, alpha_current, torque_current);
        end
    catch
        fprintf('\\n TCP transmission error at step %d\\n', k);
        break;
    end
end

%% CLEANUP
try
    fclose(tcpipClient);
    delete(tcpipClient);
catch
end
fprintf('\\n===============================================\\n');
fprintf(' Augmented Lagrange simulation complete (EoM alpha)!\\n');
fprintf(' Transmitted: θ, ω, α_EoM, τ to C++ server\\n');
fprintf('===============================================\\n\\n');

%% UPDATED LOCAL FUNCTION: Augmented Lagrange ODE (EoM version)
function dydt = augmentedLagrangeODE_EoM(y, m_A, g, L, I_theta)
    q = y(1:3); 
    qdot = y(4:6);
    
    % PRIORITY: Use simple pendulum EoM for alpha
    alpha_EoM = -(m_A * g * (L/2) * sin(q(3))) / I_theta;
    
    % For completeness, compute translational accelerations via constraints
    M = diag([m_A, m_A, I_theta]);
    Cq = [1, 0, (L/2)*sin(q(3)); 0, 1, (L/2)*cos(q(3))];
    Qe = [0; -m_A*g; 0];
    Qv_theta = -(qdot(3)^2)*(L/2)*cos(q(3)); 
    Qv = [0; 0; Qv_theta];
    Cqq_qdot2 = [qdot(3)^2*(L/2)*cos(q(3)); 
                 -qdot(3)^2*(L/2)*sin(q(3))];
    A = [M, Cq'; Cq, zeros(2)]; 
    b = [Qe + Qv; -Cqq_qdot2];
    sol = A\b; 
    
    qddot = [sol(1); sol(2); alpha_EoM]; % Override angular acc with EoM
    
    dydt = [qdot; qddot];
end


