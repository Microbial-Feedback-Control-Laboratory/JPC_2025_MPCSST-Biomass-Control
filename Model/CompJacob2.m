function FJ=CompJacob2(x0,par,u, flag_augm)
% Inputs:
% x0 initial state vect, linearization point
% par struct contains the model's parameters, can be Filter.augm.par, if so
% flag_augm is 1
% u input vector
% flag_augm is 1 if par is the updated parameters' struct

% get parameters and
kd     = par.kd;
Q = u(1);
Fin = u(2:end-1);
Fin = Fin(:);
Fout = u(end);

Sin = zeros(1,length(Fin));
Sin(1)= par.Sin;

Y_XSinv   = par.Y_XSinv;
mu_max = par.mu_max;
% if flag_augm==0, i.e. there is no augmentation, then we compute the jacobian of a normal model,
% otherwise we consider it as an extended model 
if exist('flag_augm','var') && flag_augm % flag_augm==1
     % Y_XSinv_fixed   = par.Y_XSinv_fixed;
    FJ=jacobianext(@Extended_Monod_model_in,x0);
else % flag_augm is false or not sent to the function:
     % Parameter inputs 
    mu = par.mu_max;
    Ks     = par.Ks;   % g L^-1
    Y_XCO2inv = par.Y_CO2Xinv;
    FJ=jacobianext(@dilution_Monod_model_in,x0);
end




    function dx=Extended_Monod_model_in(x)

        % States
        V = x(1); X = x(2); S = x(3); CO2 = x(4);
       
        % The inputs are Fin, Fout
        % Differential equations:
        dV   = sum(Fin)  - Fout;
        dX   = -X.*(sum(Fin)/V) + x(5)*(S ./(x(6) + S)) .*X - kd .*X; %Biomass
        dS   = (Sin-S)*(Fin/V) - x(5) .*(S ./(x(6) + S)) .*X .*Y_XSinv; % Substrate
        dCO2 = -CO2.*Q + (1 .*(S ./(x(6) + S)) .*X).*x(7); %CO2
        
        % Initial empty vector
        dx = zeros(length(x),1); % This creates an empty output column vector

        % Output:
        dx(1)= dV;
        dx(2)= dX;
        dx(3)= dS;
        dx(4)= dCO2;
        dx(5) = 0;
        dx(6) = 0;
        dx(7) = 0;
    end

    function dx=dilution_Monod_model_in(x)

        % States
        V = x(1); X = x(2); S = x(3); CO2 = x(4);
        % Manipulated variables
        % FlowVar=[Instrument.Air_flow; Instrument.Volume_in; Instrument.Volume_out];
       
        % The inputs are Fin, Fout
        % Differential equations:
        dV   = sum(Fin) - Fout;
        dX   = -X.*(sum(Fin)/V) + mu*(S ./(Ks + S)) .*X - kd .*X; %Biomass
        dS   = (Sin-S)*(Fin/V) - mu .*(S ./(Ks + S)) .*X .*Y_XSinv; % Substrate
        dCO2 = -CO2.*Q + (1 .*(S ./(Ks + S)) .*X).*Y_XCO2inv; %CO2
        
        % Initial empty vector
        dx = zeros(4,1); % This creates an empty output column vector

        % Output:
        dx(1)= dV;
        dx(2)= dX;
        dx(3)= dS;
        dx(4)= dCO2;
    end
end
