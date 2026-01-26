% Read data from CSV file
data = readtable('AL_IKF_Analytical_results.csv');

% Extract columns
time = data{:, 1};           
Theta_CPP = data{:, 4};      
Theta_IKF = data{:, 6};      
Theta_MATLAB = data{:, 2};   

Omega_CPP = data{:, 5};      
Omega_IKF = data{:, 7};      
Omega_MATLAB = data{:, 3};   

% Create figure with two subplots
figure('Position', [100, 100, 1200, 800]);

% First plot: Theta values
subplot(2, 1, 1);
plot(time, Theta_CPP, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Theta_{CPP}');
hold on;
plot(time, Theta_IKF, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'Theta_{IKF}');
plot(time, Theta_MATLAB, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Theta_{MATLAB}');
hold off;
xlabel('Time (s)', 'FontSize', 12);
ylabel('Theta (rad)', 'FontSize', 12);
title('Theta Comparison: CPP vs IKF vs MATLAB', 'FontSize', 13, 'FontWeight', 'bold');
legend('FontSize', 11, 'Location', 'best');
grid on;

% Second plot: Omega values
subplot(2, 1, 2);
plot(time, Omega_CPP, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Omega_{CPP}');
hold on;
plot(time, Omega_IKF, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'Omega_{IKF}');
plot(time, Omega_MATLAB, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Omega_{MATLAB}');
hold off;
xlabel('Time (s)', 'FontSize', 12);
ylabel('Omega (rad/s)', 'FontSize', 12);
title('Omega Comparison: CPP vs IKF vs MATLAB', 'FontSize', 13, 'FontWeight', 'bold');
legend('FontSize', 11, 'Location', 'best');
grid on;
