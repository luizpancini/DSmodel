function [E_mean,ps,E] = optimize_BL(p,INPUTS,varargin)

% Handle inputs
plotting = 0;       % TF to plot coefficients 
if ~isempty(varargin)
    plotting = varargin{1};
end

%% Loop over cases
% Pre-process
if ~isrow(INPUTS.cases), INPUTS.cases = INPUTS.cases'; end
switch INPUTS.source
    case "NASA"
        INPUTS.base_name = "NASA frame ";
    case "GU"
        INPUTS.base_name = "GU case ";
end

% Initialize error sums and vectors
n_cases = length(INPUTS.cases);
E_cn_sum = 0;           E_cn_vec = nan(n_cases,1);
E_cm_sum = 0;           E_cm_vec = nan(n_cases,1);
E_cc_sum = 0;           E_cc_vec = nan(n_cases,1);
E_sum = 0;

% Loop
j = 0;
for case_now = INPUTS.cases
    
    j = j+1;
    
    %% Initialization
    % Select number of cycles to discard
    if (INPUTS.source == "NASA" && ismember(case_now,[8023,8118,8203,9022,9101,9106])) || (INPUTS.source == "GU" && ismember(case_now,[11014191,11014381]))
        n_discard = 2;
    else
        n_discard = 1;
    end
    % Initialize
    [data,params,x0,y0,ODE_options] = initialize_model(case_now,INPUTS,n_discard,1e-3,0);
    % Override with input parameter values
    [params,ps] = override_BL_params(p,params,INPUTS.model);
    
    %% ODEs solver
    switch INPUTS.model
        case "BL"
            outputs = BL_RKF45(params.tspan,x0,y0,params,ODE_options);
    end
        
    %% Interpolate data and find normalized errors
    interp = interp_data(case_now,INPUTS,data,params,outputs,0,INPUTS.error_standard);
    % Errors
    E_cn = interp.error.cn; E_cn_vec(j) = E_cn;
    E_cm = interp.error.cm; E_cm_vec(j) = E_cm;
    E_cc = interp.error.cc; E_cc_vec(j) = E_cc;
    
    %% Plots
    if plotting
        call_plotter(INPUTS,case_now,data,params,outputs,interp,1);
    end
         
    %% Increment error sums
    E_cn_sum = E_cn_sum + E_cn;
    E_cm_sum = E_cm_sum + E_cm;
    E_cc_sum = E_cc_sum + E_cc;
    
    % Cumulative error sum
    E_sum = E_sum + (INPUTS.coefs_weights(1)*E_cn+INPUTS.coefs_weights(2)*E_cm+INPUTS.coefs_weights(3)*E_cc)/sum(INPUTS.coefs_weights);
                                         
end

% Mean error
E_mean = E_sum/numel(INPUTS.cases);

% Pack error data into structure
E = variables2struct(struct,E_mean,E_sum,E_cn_sum,E_cm_sum,E_cc_sum,E_cn_vec,E_cm_vec,E_cc_vec);
    
end