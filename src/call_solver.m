%% Pre-process
switch INPUTS.source
    case "NASA"
        INPUTS.cases = NASA;
        INPUTS.base_name = "Frame ";
    case "GU"
        INPUTS.cases = GU;
        INPUTS.base_name = "GU run ";
    case "OSU"
        INPUTS.cases = OSU;
        INPUTS.base_name = "OSU run ";
    case "other"
        INPUTS.cases = other;
        INPUTS.base_name = "Oscillatory case ";
    case "gust_test"
        INPUTS.model = INPUTS.model + "gust";
        INPUTS.cases = gust_test;
        INPUTS.base_name = "Gust test case ";
    case "flap_test"
        INPUTS.model = INPUTS.model + "flap";
        INPUTS.cases = flap_test;
        INPUTS.base_name = "Flap test case ";
    case "tvfs_test"
        INPUTS.model = INPUTS.model + "tvfs";
        INPUTS.cases = tvfs_test;
        INPUTS.base_name = "TVFS test case "; 
    case "hargen_test"
        INPUTS.model = INPUTS.model + "hargen";
        INPUTS.cases = hargen_test;
        INPUTS.base_name = "General harmonic test case ";    
    otherwise
        error("Source '" + INPUTS.source + "' not available")
end

% Loop over cases
counter = 0; 
for case_now = INPUTS.cases(:)'

% Update case counter and change data type if necessary
counter = counter+1;
if iscell(case_now), case_now = case_now{1}; end

%% Initialization 
[data,params,x0,y0,RKF_options] = initialize_model(case_now,INPUTS);

%% Solve system and get outputs
tic; 
switch INPUTS.model
    case "BL"
        outputs = BL_RKF45(params.tspan,x0,y0,params,RKF_options);        
    case "BLT"
        outputs = BLT_RKF45(params.tspan,x0,y0,params,RKF_options);         
    case {"BLS","BLSLR"}
        outputs = BLS_RKF45(params.tspan,x0,y0,params,RKF_options);
    case "BLO"
        outputs = BLO_RKF45(params.tspan,x0,y0,params,RKF_options); 
    case "BLG"
        outputs = BLG_RKF45(params.tspan,x0,y0,params,RKF_options);
    case "BLgust"
        outputs = BLgust_RKF45(params.tspan,x0,y0,params,RKF_options); 
    case "BLTgust"
        outputs = BLTgust_RKF45(params.tspan,x0,y0,params,RKF_options);
    case "BLOgust"
        outputs = BLOgust_RKF45(params.tspan,x0,y0,params,RKF_options);
    case "BLTflap"
        outputs = BLTflap_RKF45(params.tspan,x0,y0,params,RKF_options);
    case "BLtvfs"
        outputs = BLtvfs_RKF45(params.tspan,x0,y0,params,RKF_options);
    case "BLTtvfs"
        outputs = BLTtvfs_RKF45(params.tspan,x0,y0,params,RKF_options);
    case "BLhargen"
        outputs = BLhargen_RKF45(params.tspan,x0,y0,params,RKF_options); 
    case "BLThargen"
        outputs = BLThargen_RKF45(params.tspan,x0,y0,params,RKF_options); 
end
outputs.runtime = toc;
% Display runtime
fprintf(INPUTS.base_name + num2str(case_now)  + ", runtime: " + num2str(outputs.runtime) + " s");

%% Interpolate data and find normalized errors with respect to reference data
interp = interp_data(case_now,INPUTS,data,params,outputs,1,INPUTS.error_standard);

%% Plots
if INPUTS.plot_sep_figs
    % Separate figures
    h = call_plotter_sepfigs(INPUTS,case_now,data,params,outputs,interp,INPUTS.save_folder,INPUTS.figures_extension);
else
    % Tabbed figures
    h = call_plotter(INPUTS,case_now,data,params,outputs,interp);
end

%% Save all of current case's data
OUTPUTS(counter) = variables2struct(struct,data,params,outputs,interp,h);
% Prepare to print on next line
fprintf('\n');

end
% Mean errors across all cases
if isfield(OUTPUTS(1).interp,'error')
    cn_error_sum = 0; cm_error_sum = 0; cc_error_sum = 0;
    for i=1:length(OUTPUTS)
        cn_error_sum = cn_error_sum+OUTPUTS(i).interp.error.cn;
        cm_error_sum = cm_error_sum+OUTPUTS(i).interp.error.cm;
        cc_error_sum = cc_error_sum+OUTPUTS(i).interp.error.cc;
    end
    mean_cn_error = cn_error_sum/length(OUTPUTS); 
    mean_cm_error = cm_error_sum/length(OUTPUTS);
    mean_cc_error = cc_error_sum/length(OUTPUTS);
    str = "Mean errors: cn = " + num2str(mean_cn_error,'%.4f') + ", cm = " + num2str(mean_cm_error,'%.4f') + ", cc = " + num2str(mean_cc_error,'%.4f');
    disp(str);
end
% Clear workspace
clearvars -except INPUTS OUTPUTS mean_cn_NRMSE mean_cm_NRMSE mean_cc_NRMSE total_runtime