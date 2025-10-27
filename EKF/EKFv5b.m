% EKF implementation
function [Filter] = EKFv5b(Filter, meas, par, iteration, u, tspan)
% Inputs: VarMem
% VarMem struct
% Filter struct
% meas vector measurements normally [V X nan CO2]
% par struct contains the model's parameters
% iteration scalar,  current iteration
% u input vector
% tspan vector, integrating time

% Outputs:
% Updated Filter struct

if iteration == 1
    Filter.xhat = par.x0; %1x4
    Filter.P = par.P0;    
else
N = length(Filter.xhat); % number of states

% Take the previous 
x_hat0 = Filter.xhat;
P = Filter.P;
H = par.H;
R = par.R;
Q = par.Q;

if length(meas) == 4 && ~isnan(meas(3))
   R(4,4) = R(3,3);
   R(3,3) = par.R_ekf_infrequent;

   H = eye(4);
end

% Identity Matrix
I = eye(N);
    
% Model simulation
[x] = CompModel(tspan, x_hat0, u, par);

% No negative state values!
x(x(:)<0) = 0.001;
x_bar = x(end,:)';

% Measurements
z1 = meas(1);
z2 = meas(2);
z3 = meas(3);
z4 = meas(4);

if isnan(z3)
    % Measurements
    z_m = [z1; z2; z4];
else    % if HPLC sample was taken
    z_m = [z1; z2; z3; z4];

end

Pk_1 = P; % this is P(k-1|k-1)
   
FJ = CompJacob2(x_hat0,par,u,false);
% Use the code Jacobian EKF to confirm the Jacobian of the differential
% equations with respect to the states
Lk_1 = 1;
Ts = tspan(end)-tspan(1);
% Ts = par.iterationTime/3600;
%% time updates eqs
phi = eye(N) + Ts*FJ;
Pk1 = phi*Pk_1*phi'+ Lk_1 * Q * Lk_1'; % this is P(k|k-1)
%   Pk1 = FJ*Pk_1*FJ'+Q; % this is P(k|k-1)

% Compute update:
yk = z_m-H*x_bar;
Sk = H *Pk1*H' +R;
Kk = Pk1*H' *inv(Sk);

%K = P * H' * inv(H * P * H' + R);

% x_hat2 is x_hat(k|k), the corrected
xhat2 = x_bar + Kk * yk;

% not possible negative states
xhat2(xhat2<0) = 0.0001;


% Replace NaN with a value close to zero
xhat2(isnan(xhat2)) = 0.0001;
Pkk = (I - Kk * H) * Pk1 * (I - Kk * H)' + Kk * R * Kk';
Filter.xhat = xhat2;
Filter.P = Pkk;
end
end