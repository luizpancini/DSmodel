function [E_sum,E,ps] = optimize_BL_plunge(p,INPUTS,varargin)

% Handle inputs
plotting = 0;           % TF to plot  
if ~isempty(varargin)
    plotting = varargin{1};
end

%% Loop over cases
% Pre-process
if ~isrow(INPUTS.cases), INPUTS.cases = INPUTS.cases'; end
INPUTS.base_name = "Oscillation test case ";

% Initialize error sum and error arrays
E_sum = 0; 
E_cn_vec = nan(length(INPUTS.cases),1);
E_cm_vec = nan(length(INPUTS.cases),1);

% Loop
j = 0;
for case_now = INPUTS.cases
    
    j = j+1;
    
    %% Initialization
    [data,params,x0,y0,ODE_options] = initialize_model(case_now,INPUTS,1,1e-3,0);
    % Set convection gust parameters
    [params,ps] = override_BL_params(p,params,INPUTS.model);
    
    %% ODEs solver
    outputs = BL_RKF45(params.tspan,x0,y0,params,ODE_options);
        
    %% Interpolate data and find normalized errors
    interp = interp_data(case_now,INPUTS,data,params,outputs,0);
    % Errors in aerodynamic coefficients
    E_cn = interp.NRMSE.cn; E_cn_vec(j) = E_cn;
    E_cm = interp.NRMSE.cm; E_cm_vec(j) = E_cm;
    % Normalized cumulative error sum
    E_sum = E_sum + INPUTS.cases_weights(j)*(INPUTS.coefs_weights.cn*E_cn+INPUTS.coefs_weights.cm*E_cm)/(INPUTS.coefs_weights.cn+INPUTS.coefs_weights.cm);
                                  
    %% Plots
    if plotting
        call_plotter(INPUTS,case_now,data,params,outputs,interp,0);
    end 
    
end

% Pack error data into structure
E = variables2struct(struct,E_sum,E_cn_vec,E_cm_vec);
       
end