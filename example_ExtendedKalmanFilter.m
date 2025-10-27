%% Extended Kalman Filter
% %% Author: Lipe Carmel and Giacomo Sartori
% %% Affiliation: Department of Chemical Engineering, Norwegian University of Science and Technolog
% %% Email: lipe.carmel@ntnu.no
% %% Date: October 24, 2025
%
% %% Description:
% %   This script is a brief tutorial on the use of the EKFv5b function.

clc; clearvars; close all;
% Add folders and subfolders to path
current_dir = fileparts(mfilename('fullpath'));
addpath(genpath(current_dir)) 

%% Model definition
par.Sin = 70;
par.mu_max = 0.19445;
par.mu = par.mu_max;
par.Ks = 0.007;
par.kd = 0.006;
par.Y_XSinv = 1/0.42042;
par.Y_CO2X = 0.54308;
par.Y_CO2Xinv = 1/0.54308;
Q = 2;
state_dot = @(x, u) dilution_Monod_model(0, x, [Q, u(:)'], par);

%% Tuning parameters used in the paper
par.P0 = diag([0.01 0.037 0.0297 0.01]); % Initial posteriori covariance matrix
par.Q = diag([0.01 0.01 0.05 0.0015]);   % Process noise covariance matrix 
par.R = diag([0.001 (0.01) (0.05)]);     % Measurement noise covariance matrix 

par.R_ekf_infrequent = 0.1; % Used only when substrate measurements are available

%% Initializing
Filter = struct;
% The initial state vector will be used as the first xhat
par.x0 = [1 10 5 0.004];

% The measurement matrix is initialized assuming no substrate measurements.
% The estimator function EKFv5b() changes H automatically to handle
% substrate measurements when available.
par.H = [ 1     0     0     0
          0     1     0     0
          0     0     0     1 ];


% On the first iteration the filter is initialized and xhat is set to the
% initial state in par
iteration = 1;
[Filter] = EKFv5b(Filter, [], par, iteration, [], []);
x_estimate = Filter.xhat; 
disp(x_estimate)

%% Measurement without substrate
% On the following iterations the current measurement and the previous
% input (assumed piecewise constant) are also required
iteration = 2;

% The measurement vector can contain a NaN in place of the substrate value
measurement = [0.9 10 NaN 0.004];

% tspan is used in the forward Euler method. It should be small in the time scale of the system
tspan = [0 1/60]; 

% Input vector
u = [0, 0, 0.1];

[Filter] = EKFv5b(Filter, measurement, par, iteration, u, tspan);
x_estimate = Filter.xhat; 
disp(x_estimate)

%% Measurement with substrate
iteration = 3;

% If a substrate measurement is available, use it instead of NaN
% The measurement matrix remains unchanged in par
measurement = [0.9 10 10 0.004];

[Filter] = EKFv5b(Filter, measurement, par, iteration, u, tspan);
x_estimate = Filter.xhat; 
disp(x_estimate)