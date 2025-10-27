function y_sol = CompModel(tspan, xk, uk, param,adaptPar,different_model, different_solver)
    % CompModel Simulates a dynamic system model over a specified time span.
    %
    % Syntax:
    %   y_sol = CompModel(tspan, xk, uk, param)
    %   y_sol = CompModel(tspan, xk, uk, param, different_model)
    %   y_sol = CompModel(tspan, xk, uk, param, different_model, different_solver)
    %
    % Description:
    %   CompModel integrates a system of ordinary differential equations (ODEs)
    %   using a specified model and solver. It allows for flexibility in selecting
    %   different models and numerical solvers as needed.
    %
    % Inputs:
    %   tspan            - [1x2 vector] Time span for the simulation.
    %
    %   xk               - [Nx1 vector] Initial state vector of the system at the start time t0.
    %
    %   uk               - [Vector] Control input(s) applied to the system.
    %                      - Assumed to be constant over the entire simulation period.
    %
    %   param            - [Struct] All other parameters required by the model.
    %
    %   adaptPar         -  [Struct] adaptive parameters VarMem.mu and
    %                       VarMem.beta. This should recompute par.mu and
    %                       par.YCO2inv according to the adaptive model. 
    %                       if you send adaptpar, it must contain all the
    %                       necessary parameters for the adaptive models
    %
    %   different_model  - (Optional) [Function Handle] Alternative model function to use instead of the default.
    %                      - Must be in the form dydt = altModel(t, x, uk, param).
    %                      - SHOULD ONLY BE USED IN EXCEPTIONAL CASES.
    %
    %   different_solver - (Optional) [Function Handle] Alternative ODE solver to use instead of the default.
    %                      - Can be any MATLAB ODE solver like @ode45, @ode23s, etc.
    %                      - SHOULD ONLY BE USED IN EXCEPTIONAL CASES.
    %
    % Outputs:
    %   y_sol            - [1xM vector] Final state vector of the system at the end of the simulation period.
    %

if nargin<5
    adaptPar=0;
end
    % Model selection
    % Assumes ode model with constant input over
   % model = @(t,x) temp_model(t, x, uk, param);
    model = @(t,x) dilution_Monod_model(t, x, uk, param,adaptPar);
   
    % model = @(t,x) OTHER OPTION(t,x,u,param);
    % model = @(t,x) ANOTHER OPTION(t,x,u,param);

    if nargin > 5
        if isa(different_model,'function_handle')
            % If another model was provided, then overwrite
            % the default.
            % This is premitted to allow the function to be used for
            % comparing different alternatives and running then in
            % parallel. For example: Monod + New mechanistic model + NN
           
            model = different_model;
        end
    end


    % Solver selection
    % solver = @ode15s;
    solver = @RKF45_book;
    % model = @(t,x) OTHER OPTION(t,x,u,param);
    % model = @(t,x) ANOTHER OPTION(t,x,u,param);

    if nargin > 6
        if isa(different_solver,'function_handle')
            % If another solver was provided, for example, to run this
            % function multiple times while comparing solvers then overwrite
            % the default.
            solver = different_solver;
        end
    end

    [~, y] = solver(@(t, x) model(t, x), tspan, xk);
    y_sol =  y(end,:);

end