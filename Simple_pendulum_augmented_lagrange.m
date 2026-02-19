%% PC2_Augmented_Lagrange_EoM - WITH SENSOR NOISE (changed the notation to suraj´s paper notations)
clear; clc; close all;

%% SYSTEM PARAMETERS (MATCH C++ SERVER) - Matches paper inertial parameters
L = 1.0;        % Length of body [m]
w = 0.1;        % Width [m]
h = 0.1;        % Height [m]
rho = 7850;     % Density [kg/m^3]
g = 8;          % Gravity [m/s^2] - Matches C++

% Mass and inertia (paper: m, J_θ)
m = rho * L * w * h;      
I_theta = m * L^2 / 12;   

%% SENSOR NOISE PARAMETERS 
sigma_z     = 1e-3;     % Position noise std dev [rad] ~0.057°
sigma_zdot  = 1e-2;     % Velocity noise std dev [rad/s]
sigma_zddot = 0.05;     % Acceleration noise std dev [rad/s²]

% Generate noise sequences (white Gaussian noise)
rng(42);
noise_z     = sigma_z     * randn(1, 10000);
noise_zdot  = sigma_zdot  * randn(1, 10000);
noise_zddot = sigma_zddot * randn(1, 10000);

fprintf('===============================================\n');
fprintf('AUGMENTED LAGRANGE - SIMPLE PENDULUM (SENSOR NOISE)\n');
fprintf('===============================================\n');
fprintf('Length L = %.3f m\n', L);
fprintf('Mass m = %.3f kg\n', m);
fprintf('Inertia I_θ = %.6f kg·m²\n', I_theta);
fprintf('Gravity g = %.2f m/s²\n', g);
fprintf('\nSENSOR NOISE:\n');
fprintf('  sigma_z_tilde     = %.1e rad (%.3f°)\n', sigma_z, sigma_z*180/pi);
fprintf('  sigma_zdot_tilde  = %.1e rad/s\n', sigma_zdot);
fprintf('  sigma_zddot_tilde = %.3f rad/s²\n', sigma_zddot);
fprintf('===============================================\n\n');

%% SIMULATION PARAMETERS
dt = 0.001;     % Timestep [s] (1 kHz)
t_end = 10.0;   % Duration [s]
t = (0:dt:t_end)';
N = length(t);

%% INITIAL CONDITIONS 
q0 = [0.1; 0.0; pi/6];      
qdot0 = [0; 0; 0.5];        
y0 = [q0; qdot0];           
fprintf('Initial Conditions (Paper: z(0), \dot z(0)):\n');
fprintf(' q(0) = [%.3f, %.3f, %.4f rad]\n', q0);
fprintf(' qdot(0) = [%.3f, %.3f, %.3f rad/s]\n', qdot0);
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
y = zeros(6, N);            
y_noisy = zeros(6, N);      % Store noisy measurements
y(:,1) = y0;
y_noisy(:,1) = y0;
noise_idx = 1;

%% MAIN AUGMENTED LAGRANGE SIMULATION LOOP WITH SENSOR NOISE

tic;
fprintf('===============================================\n');
fprintf('AUGMENTED LAGRANGE SIMULATION STARTED (w/ NOISE)\n');
fprintf('===============================================\n');
for k = 1:N-1
    % Current state (Paper: q, \dot q)
    q = y(1:3,k);
    qdot = y(4:6,k);
    
    % ===============================================
    % STEP 1: ASSEMBLE MASS MATRIX M (3x3)
    % ===============================================
    M = diag([m, m, I_theta]);
    
    % ===============================================
    % STEP 2: CONSTRAINT JACOBIAN Φ_q (2x3) 
    % ===============================================
    Phi_q = [1, 0, (L/2)*sin(q(3));
             0, 1, (L/2)*cos(q(3))];
    
    % ===============================================
    % STEP 3: EXTERNAL FORCES Q (3x1) - Paper: Q
    % ===============================================
    Q = [0; -m*g; 0];
    
    % ===============================================
    % STEP 4: QUADRATIC VELOCITY Q_v (3x1)
    % ===============================================
    Qv_theta = -(qdot(3)^2)*(L/2)*cos(q(3));
    Q_v = [0; 0; Qv_theta];
    
    % ===============================================
    % STEP 5: CONSTRAINT ACCELERATION 
    % ===============================================
    Phi_qq_qdot2 = [qdot(3)^2*(L/2)*cos(q(3));
                    -qdot(3)^2*(L/2)*sin(q(3))];
    
    % ===============================================
    % STEP 6: SOLVE AUGMENTED SYSTEM 
    % ===============================================
    A_aug = [M, Phi_q';
             Phi_q, zeros(2)];
    b_aug = [Q + Q_v; -Phi_qq_qdot2];
    sol = A_aug \ b_aug;
    
    qddot = sol(1:3);       
    lambda = sol(4:5);      
    
    zddot = qddot(3);        
    
    % ===============================================
    % STEP 7: RK4 INTEGRATION (TRUE DYNAMICS)
    % ===============================================
    k1 = augmentedLagrangeODE(y(:,k), m, g, L, I_theta);
    k2 = augmentedLagrangeODE(y(:,k) + dt/2*k1, m, g, L, I_theta);
    k3 = augmentedLagrangeODE(y(:,k) + dt/2*k2, m, g, L, I_theta);
    k4 = augmentedLagrangeODE(y(:,k) + dt*k3, m, g, L, I_theta);
    
    y_next = y(:,k) + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    y(:,k+1) = y_next;
    
    % ===============================================
    % STEP 8: ADD SENSOR NOISE TO MEASUREMENTS 
    % ===============================================
    z_true = y_next(3);     % True joint position: z
    zdot_true = y_next(6);  % True joint velocity: \dot z
    zddot_true = zddot;     % True joint acceleration: \ddot z
    
    % Apply sensor noise (Paper: z, dot z, ddot z)
    z_tilde = z_true + noise_z(noise_idx);
    zdot_tilde = zdot_true + noise_zdot(noise_idx);
    zddot_tilde = zddot_true + noise_zddot(noise_idx);
    
    % Store noisy measurements
    y_noisy(3,k+1) = z_tilde;
    y_noisy(6,k+1) = zdot_tilde;
    
    noise_idx = noise_idx + 1;
    if noise_idx > length(noise_z)
        noise_idx = 1; % Wrap around
    end
    
    % ===============================================
    % STEP 9: TCP TRANSMISSION 
    % ===============================================
    try
        fwrite(tcpipClient, z_tilde, 'double');    % Noisy position \tilde z
        fwrite(tcpipClient, zdot_tilde, 'double'); % Noisy velocity \tilde{\dot z}
        fwrite(tcpipClient, zddot_tilde, 'double');% Noisy acceleration \tilde{\ddot z}
        
        if mod(k, 1000) == 0 % Progress every 1s
            fprintf('Step %d: z=%.4f/%.4f, zdot=%.4f/%.4f, zddot=%.4f/%.4f\n', ...
                k, z_tilde, z_true, zdot_tilde, zdot_true, zddot_tilde, zddot_true);
        end
    catch
        fprintf('\n✗ TCP transmission error at step %d\n', k);
        break;
    end
end

%% CLEANUP
try
    fclose(tcpipClient);
    delete(tcpipClient);
catch
end
elapsed = toc;
fprintf('\n===============================================\n');
fprintf('  Simulation complete with sensor noise!\n');
fprintf('  Duration: %.2f s (realtime factor: %.2fx)\n', elapsed, t_end/elapsed);
fprintf('  Transmitted noisy: z_tilde, zdot_tilde, zddot_tilde to C++ server\n');
fprintf('===============================================\n\n');

%% PLOT RESULTS WITH NOISE COMPARISON
figure('Name', 'MATLAB Results: True vs Noisy Measurements', 'Position', [100 100 1400 900]);

% True vs Noisy states
subplot(3,1,1);
plot(t, y(3,:), 'b-', 'LineWidth', 2, 'DisplayName', 'True z');
hold on;
plot(t, y_noisy(3,:), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Noisy z_tilde');
grid on; legend('Location', 'best');
xlabel('Time [s]'); ylabel('z [rad]');
title('Pendulum Angle: True vs Sensor Measurements');

subplot(3,1,2);
plot(t, y(6,:), 'b-', 'LineWidth', 2, 'DisplayName', 'True zdot');
hold on;
plot(t, y_noisy(6,:), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Noisy zdot_tilde');
grid on; legend('Location', 'best');
xlabel('Time [s]'); ylabel('zdot [rad/s]');
title('Angular Velocity: True vs Sensor Measurements');

%% NOISE STATISTICS 
fprintf('NOISE STATISTICS (over %d samples):\n', N-1);
fprintf('  z noise:     RMS=%.3e rad,  max=%.3e rad\n', ...
        rms(noise_z(1:N-1)), max(abs(noise_z(1:N-1))));
fprintf('  zdot noise:  RMS=%.3e rad/s, max=%.3e rad/s\n', ...
        rms(noise_zdot(1:N-1)), max(abs(noise_zdot(1:N-1))));
fprintf('  zddot noise: RMS=%.3e rad/s², max=%.3e rad/s²\n', ...
        rms(noise_zddot(1:N-1)), max(abs(noise_zddot(1:N-1))));

%% LOCAL FUNCTION: Augmented Lagrange ODE 
function dydt = augmentedLagrangeODE(y, m, g, L, I_theta)
    q = y(1:3); 
    qdot = y(4:6);
    
    M = diag([m, m, I_theta]);
    Phi_q = [1, 0, (L/2)*sin(q(3)); 0, 1, (L/2)*cos(q(3))];
    Q = [0; -m*g; 0];
    Qv_theta = -(qdot(3)^2)*(L/2)*cos(q(3)); 
    Q_v = [0; 0; Qv_theta];
    Phi_qq_qdot2 = [qdot(3)^2*(L/2)*cos(q(3)); -qdot(3)^2*(L/2)*sin(q(3))];
    
    A_aug = [M, Phi_q'; Phi_q, zeros(2)]; 
    b_aug = [Q + Q_v; -Phi_qq_qdot2];
    sol = A_aug\b_aug; 
    
    qddot = sol(1:3);
    dydt = [qdot; qddot];
end
