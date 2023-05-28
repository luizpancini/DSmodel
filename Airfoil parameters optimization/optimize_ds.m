function E_sum = optimize_ds(p)

%% Inputs
OSU = [372 373 374 384 385 386 396 408 409 410 420 421 422];           
% Select input option as "NASA", "GU", "OSU"
INPUTS.source = "OSU";
% Select model (BL or BLS, for Sheng et al.'s model)
INPUTS.model = "BL";
% Gust indicial response model
INPUTS.gust_ind = "";

%% Loop over cases
% Pre-process
switch INPUTS.source
    case "NASA"
        INPUTS.cases = frame;
    case "GU"
        INPUTS.cases = GUD;
    case "OSU"
        INPUTS.cases = OSU;
end
if ~isrow(INPUTS.cases), INPUTS.cases = INPUTS.cases'; end

% Initialize error sums
E_cn = 0;
E_cm = 0;
E_cc = 0;
E_sum = 0;

% Loop
j = 0;
for case_now = INPUTS.cases
    
    j = j+1;
    
    %% Initialization
    [data,params,x0,y0,ODE_options] = initialize_model(case_now,INPUTS);
    % Override with input parameter values
    params.alpha_ss = p(1);
    params.alpha_ds0 = p(2);
    params.r0 = p(3);
    params.Ta = p(4);
    
    %% ODEs solver
%     tic;
    switch INPUTS.model
        case "BL"
            [t,x,y,xdot] = BL_RKF45(params.tspan,x0,y0,params,ODE_options);
    end
%     run_time = toc;
%     str = input_name + num2str(case_now) + ", runtime: " + num2str(run_time) + " s"; fprintf(str);

    %% Get output variables
    switch INPUTS.model
        case "BL"
            outputs = BL_output_vars(x,y,t);
    end
    
    %% Interpolate data and find normalized errors
    interp = interp_data(INPUTS,data,params,outputs,0);
    
    %% Increment error sums
    E_cn = E_cn + interp.NRMSE.cn;
    E_cm = E_cm + interp.NRMSE.cm;
    E_cc = E_cc + interp.NRMSE.cc;
    
    E_sum = E_sum + E_cn + E_cc;
end

end