%% PC2_Augmented_Lagrange_Clean_NoPlots.m
% PURE AUGMENTED LAGRANGE METHOD IMPLEMENTATION
% NO PLOTTING - Clean version for thesis

clear; clc; close all;

%% SYSTEM PARAMETERS 
L = 1.0;           % Length of body [m]
w = 0.1;           % Width [m]
h = 0.1;           % Height [m]
rho = 7850;        % Density [kg/m^3]
g = 9.81;          % Gravity [m/s^2]

% Mass and inertia
m_A = rho * L * w * h;                    % Mass m_A
I_theta = m_A * L^2 / 12;                 % Moment of inertia about COM

fprintf('===============================================\n');
fprintf('AUGMENTED LAGRANGE - SIMPLE PENDULUM\n');
fprintf('===============================================\n');
fprintf('Length l     = %.3f m\n', L);
fprintf('Mass m_A     = %.3f kg\n', m_A);
fprintf('Inertia I_θ  = %.6f kg·m²\n', I_theta);
fprintf('===============================================\n\n');

%% SIMULATION PARAMETERS
dt = 0.001;              % Timestep [s] (1 kHz)
t_end = 10.0;            % Duration [s]
t = (0:dt:t_end)';       % Time vector
N = length(t);

%% INITIAL CONDITIONS
% q = [R_x, R_y, theta], qdot = [R_x_dot, R_y_dot, theta_dot]
q0 = [0.1; 0.0; pi/6];   % Initial: small x offset, theta=30°
qdot0 = [0; 0; 0.5];     % Initial velocities
y0 = [q0; qdot0];        % Full state [6x1]

fprintf('Initial Conditions:\n');
fprintf(' q(0) = [%.3f, %.3f, %.3f] rad\n', q0);
fprintf(' qdot(0) = [%.3f, %.3f, %.3f] rad/s\n', qdot0);
fprintf('\n');

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
    fprintf('✓ Connected to C++ Augmented Lagrange server!\n\n');
catch ME
    fprintf('✗ Connection failed: %s\n', ME.message);
    return;
end

%% SIMULATION STORAGE
y = zeros(6, N);     % [R_x, R_y, θ, R_x_dot, R_y_dot, θ_dot]
y(:,1) = y0;
lambda_log = zeros(2, N);  % Lagrange multipliers

%% MAIN AUGMENTED LAGRANGE SIMULATION LOOP
tic;
fprintf('===============================================\n');
fprintf('AUGMENTED LAGRANGE SIMULATION STARTED\n');
fprintf('===============================================\n');
fprintf('Time | θ [rad] | ω [rad/s] | λ₁ [N] | λ₂ [N]\n');
fprintf('------------------------------------------------\n');

for k = 1:N-1
    
    % Current state
    q = y(1:3,k);
    qdot = y(4:6,k);
    
    % ===============================================
    % STEP 1: ASSEMBLE MASS MATRIX M (3x3)
    % ===============================================
    M = diag([m_A, m_A, I_theta]);
    
    % ===============================================
    % STEP 2: CONSTRAINT JACOBIAN C_q (2x3)
    % ===============================================
    Cq = [1,              0,              (L/2)*sin(q(3));
          0,              1,              (L/2)*cos(q(3))];
    
    % ===============================================
    % STEP 3: GENERALIZED FORCES Q_e (3x1)
    % ===============================================
    Qe = [0; -m_A*g; 0];
    
    % ===============================================
    % STEP 4: QUADRATIC VELOCITY Q_v (3x1)
    % ===============================================
    % From C_q * qdot * qdot term
    Qv_theta = -(qdot(3)^2)*(L/2)*cos(q(3));
    Qv = [0; 0; Qv_theta];
    
    % ===============================================
    % STEP 5: CONSTRAINT ACCELERATION TERMS
    % ===============================================
    % C_qq * qdot * qdot (2x1)
    Cqq_qdot2 = [qdot(3)^2*(L/2)*cos(q(3));
                 -qdot(3)^2*(L/2)*sin(q(3))];
    
    % ===============================================
    % STEP 6: SOLVE 5x5 AUGMENTED SYSTEM
    % ===============================================
    % [M     Cq^T] [qddot]   [Qe + Qv]
    % [Cq     0  ] [lambda] = [-Cqq_qdot2]
    
    A = [M,               Cq';
         Cq,              zeros(2)];
    
    b = [Qe + Qv; -Cqq_qdot2];
    
    sol = A \ b;  % Solve 5x1 system
    
    qddot = sol(1:3);
    lambda = sol(4:5);  % Lagrange multipliers (constraint forces)
    
    % Store Lagrange multipliers
    lambda_log(:,k) = lambda;
    
    % ===============================================
    % STEP 7: RK4 INTEGRATION (full 6-state)
    % ===============================================
    k1 = augmentedLagrangeODE(y(:,k), m_A, g, L, I_theta);
    k2 = augmentedLagrangeODE(y(:,k) + dt/2*k1, m_A, g, L, I_theta);
    k3 = augmentedLagrangeODE(y(:,k) + dt/2*k2, m_A, g, L, I_theta);
    k4 = augmentedLagrangeODE(y(:,k) + dt*k3, m_A, g, L, I_theta);
    
    y_next = y(:,k) + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    y(:,k+1) = y_next;
    
    % ===============================================
    % STEP 8: TCP TRANSMISSION (send angle & angular velocity)
    % ===============================================
    theta_current = y_next(3);
    omega_current = y_next(6);
    
    try
        fwrite(tcpipClient, theta_current, 'double');
        fwrite(tcpipClient, omega_current, 'double');
    catch
        fprintf('\n✗ TCP transmission error\n');
        break;
    end
    
    % ===============================================
    % STEP 9: DISPLAY (every 100 samples = 100 ms)
    % ===============================================
    if mod(k, 100) == 0
        t_current = t(k);
        
        % Console display
        fprintf('%.3f | %+.6f | %+.6f | %+.2f | %+.2f\n', ...
            t_current, theta_current, omega_current, lambda);
    end
    
    % Real-time pacing
    while toc < t(k+1)
        pause(0.00001);
    end
end

%% CLEANUP
try
    fclose(tcpipClient);
    delete(tcpipClient);
catch
end

fprintf('\n===============================================\n');
fprintf(' Augmented Lagrange simulation complete!\n');
fprintf(' All measurements transmitted to C++ server\n');
fprintf('===============================================\n\n');

%% =========================================================
% LOCAL FUNCTION: Augmented Lagrange ODE
% =========================================================
function dydt = augmentedLagrangeODE(y, m_A, g, L, I_theta)
    % y = [R_x, R_y, theta, R_x_dot, R_y_dot, theta_dot]
    
    q = y(1:3);
    qdot = y(4:6);
    
    % Mass matrix M (3x3)
    M = diag([m_A, m_A, I_theta]);
    
    % Constraint Jacobian C_q (2x3)
    Cq = [1,              0,              (L/2)*sin(q(3));
          0,              1,              (L/2)*cos(q(3))];
    
    % Generalized forces Q_e
    Qe = [0; -m_A*g; 0];
    
    % Quadratic velocity Q_v
    Qv_theta = -(qdot(3)^2)*(L/2)*cos(q(3));
    Qv = [0; 0; Qv_theta];
    
    % Constraint acceleration terms
    Cqq_qdot2 = [qdot(3)^2*(L/2)*cos(q(3));
                 -qdot(3)^2*(L/2)*sin(q(3))];
    
    % Solve augmented system
    A = [M, Cq'; Cq, zeros(2)];
    b = [Qe + Qv; -Cqq_qdot2];
    sol = A\b;
    
    qddot = sol(1:3);
    
    % Return derivatives
    dydt = zeros(6,1);
    dydt(1:3) = qdot;      % Position derivatives (dq/dt)
    dydt(4:6) = qddot;     % Velocity derivatives (dqdot/dt)
end