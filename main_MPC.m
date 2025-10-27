clc; clear all; close all;

% Add folders and subfolders to path
current_dir = fileparts(mfilename('fullpath'));
addpath(genpath(current_dir)) 

%% System sampling rate
dt = 1/60;

%% Model definition
par.Sin = 70;
par.mu_max = 0.19445;
par.mu = par.mu_max;
par.Ks = 0.007;
par.kd = 0.006;
par.Y_XSinv = 1/0.42042;
par.Y_CO2X = 0.54308;
par.Y_CO2Xinv = 1/0.54308;

state_dot = @(x, u) dilution_Monod_model(0, x, [2, u(:)'], par);

%% MPC definition

% MPC_mono_dilution object requires a 4 x 3 process model to be constructed
% Default setpoints, tuning parameters, and constraints are predefined.
% Control actions are subjected to non-negativity constraints.
MPC = MPC_mono_dilution(state_dot);

% Properties can be overwritten
MPC.Ts = dt; % Sampling time in hours

% Properties that modify the constraint matrices of the optimization
% require the constraints() method to be enforced
MPC.p = 60;
MPC.m = 6;
MPC.constraints();

% Targets are set as properties
MPC.Xsp = 10;
MPC.Vsp = 1;

% Parallel computation can be faster on most machines
% if not(MPC.optimizer_options.UseParallel)
%     % setup parallel
%     delete(gcp('nocreate'))
%     c = parcluster('Processes');
%     parpool(2, 'IdleTimeout', inf);
% 
%     MPC.optimizer_options.UseParallel = true;
%     MPC.optimizer_options.MaxIterations = 100;
% end

% To calculate control actions use MPC.solve(x, u) where x and u are row
% vectors representing the current state and previous action.
% See:
help MPC.solve

%% Creating the simulated process
plant = @(x, u) dilution_Monod_model(0, x, [2, u(:)'], par);

%% Initial conditions

% Duration of the simulation
tf = 10;                % h

V_0 = 1;
X_0 = 5;
S_0 = 5;
CO2_0 = 0.0049;
x0_plant = [V_0, X_0, S_0, CO2_0];

uk = [0, 0, 0];

nu = length(uk);
nx = length(x0_plant);


%% EKF initialization
par.P0 = diag([0.01 0.037 0.0297 0.01]);
par.Q = diag([0.01 0.01 0.05 0.0015]); 
par.R = diag([0.001 (0.01) (0.05)]);
par.R_ekf_infrequent = 0.1;

par.x0 =  x0_plant;
Filter.xhat = x0_plant;
Filter.P = par.P0;

par.H = [ 1     0     0     0
          0     1     0     0
          0     0     0     1 ];

%% Simulation setup
tspan= [0 dt]; % from t to t+dt, or just 0 to dt because the ODE is time invariant

N = ceil(tf/dt) + 1;   % number of time samples + initial condition

% Initialization
U = zeros(N, nu);
Y_plant = zeros(N, nx);
Y_measurements = Y_plant;
Y_estimate= Y_plant;
Y_model = Y_plant;
Y_sp = zeros(N, nx - 1); % CO2 is open-loop


U(1, :) = uk;


opts = odeset('NonNegative', [2 3]);

%% Simulation
timer = tic;
for i = 1 : N
    fprintf('Simulated: %.1f %% \n', i/N*100)
    fprintf('Time elapsed: %.1f minutes \n', toc(timer)/60)

    Y_plant(i, :) = x0_plant;

    %% Measure the plant
    if i == 1
        x0_measured = x0_plant; % perfect initial measurement
    else
        x0_measured = measure(x0_plant);
        if mod(i, floor(0.5/dt)) == 0
            x0_measured(3) = x0_plant(3); % substrate measurement
        end
    end
    Y_measurements(i,:) = x0_measured;

    %% Estimate the states (uk is uk-1)
    [Filter] = EKFv5b(Filter, x0_measured, par, i, uk, tspan);
    x0_estimate = Filter.xhat;
    Y_estimate(i,:) = x0_estimate;

    %% Calculate control action
    Kp = 2;
    Ssp = Kp*(MPC.Xsp - x0_estimate(2));
    Ssp = min([3 Ssp]);
    MPC.Ssp = Ssp;

    if x0_estimate(2) > MPC.Xsp
        MPC.dilution = 0; % dilution is allowed
    else
        % Toggle slacked contraint
        MPC.dilution = 1e7;
    end

    Y_sp(i,:) = [MPC.Vsp, MPC.Xsp, MPC.Ssp];
    if mod(i, 2) ==0
    uk = MPC.solve(x0_estimate(:)', uk(:)');
    end

    U(i, :) = uk;

    fprintf('\nControl action: \n')
    disp(uk)

    %% Apply control action to the process and obtain y k+1

    % Plant
    [t,y] = ode45(@(t,x) plant(x, uk), tspan, x0_plant, opts);

    %% ---------------------------- k+1 -----------------------------------
    x0_plant = y(end, :);

end

%% Results
t = 0 : dt : (i-1)*dt;

figure(1);
clf

% --- Volume ---
subplot(3,1,1);
plot(t, Y_plant(1:i,1), 'b-', 'LineWidth', 3, 'DisplayName', 'Plant'); hold on;
plot(t, Y_estimate(1:i,1), 'cyan--', 'LineWidth', 3, 'DisplayName', 'Estimate');
plot(t, Y_sp(1:i,1), 'r--', 'LineWidth', 3, 'DisplayName', 'Setpoint');
grid on; box on;
xlabel('Time (h)');
ylabel('V (L)');
legend('Location','best');
hold off;

% --- Biomass concentration ---
subplot(3,1,2);
plot(t, Y_plant(1:i,2), 'b-', 'LineWidth', 3, 'DisplayName', 'Plant'); hold on;
plot(t, Y_estimate(1:i,2), 'cyan--', 'LineWidth', 3, 'DisplayName', 'Estimate');
plot(t, Y_sp(1:i,2), 'r--', 'LineWidth', 3, 'DisplayName', 'Setpoint');
grid on; box on;
xlabel('Time (h)');
ylabel('X (g/L)');
legend('Location','best');
hold off;

% --- Substrate concentration ---
subplot(3,1,3);
plot(t, Y_plant(1:i,3), 'k-', 'LineWidth', 3, 'DisplayName', 'Plant'); hold on;
plot(t, Y_estimate(1:i,3), 'cyan--', 'LineWidth', 3, 'DisplayName', 'Estimate');
plot(t, Y_sp(1:i,3), 'r--', 'LineWidth', 3, 'DisplayName', 'Setpoint');
grid on; box on;
xlabel('Time (h)');
ylabel('S (g/L)');
legend('Location','best');
hold off;


ax = findall(gcf, 'type', 'axes');

for i = 1:length(ax)
    ax(i).FontSize = 15;
    ax(i).XLabel.FontSize = 15;
    ax(i).YLabel.FontSize = 15;
end

function x = measure(x)
    %x(1) = x(1) + normrnd(0, 0.001);
    x(2) = x(2) + normrnd(0, 0.001);
    %x(4) = x(4) + normrnd(0, 0.01);

    x(3) = nan; % Substrate is less frequently sampled

    for i = 1 : length(x)
        if x(i) < 0
            x(i) = 0;
        end
    end
end
