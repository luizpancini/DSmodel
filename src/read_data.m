function [data,params] = read_data(INPUTS,case_now,n_discard)

% Unpack
source = INPUTS.source;
model = INPUTS.model;
gust_ind = INPUTS.gust_ind;

%% Load experimental/model data and test conditions from source
switch source
    case "NASA"
        [params,data] = load_NASA(case_now);
    case "GU"
        [params,data] = load_GU(case_now);    
	case "OSU"
        [params,data] = load_OSU(case_now);
	case "other"
        [params,data] = load_other(case_now);
    case "gust_test"
        [params,data] = load_gust_test(case_now);
    case "flap_test"
        [params,data] = load_flap_test(case_now);
    case "tvfs_test"
        [params,data] = load_tvfs_test(case_now);  
    case "hargen_test"
        [params,data] = load_hargen_test(case_now);  
    otherwise
        error("Unknown data source: " + source);
end

% Unpack
airfoil = params.airfoil;
b = params.b;
U = params.U;
M = params.M;
k = params.k;

% Simulation time variables for oscillatory cases
if ~(contains(model,'gust') || (isfield(data,'source') && data.source == "gust_test"))                
    if contains(model,'flap')
        % Smallest non-zero reduced frequency
        ks = [k, params.k_f]; ks = min(ks(ks>0));
        % Time of one cycle
        t_cycle = 2*pi*b/(U*ks);
    elseif contains(model,'tvfs')
        % Smallest non-zero reduced frequency
        ks = [k, params.k_U]; ks = min(ks(ks>0));
        % Time of one cycle
        t_cycle = 2*pi*b/(U*ks);
    elseif contains(model,'gen')
        % Smallest non-zero reduced frequency
        ks = [k, params.k_h, params.k_U, params.k_f]; ks = min(ks(ks>0));
        % Time of one cycle
        if ~isempty(ks)
            t_cycle = 2*pi*b/(U*ks);
        else
            t_cycle = params.tspan(end);
            ks = 0;
        end
    else
        % Smallest non-zero reduced frequency
        ks = [k, params.k_h]; ks = min(ks(ks>0));
        % Time of one cycle
        t_cycle = 2*pi*b/(U*ks);
    end
    % Number of cycles to discard for plots (only one if conditions are quasi-steady, i.e., k < 0.01)
    is_quasi_steady = ks<0.01;
    params.n_discard = n_discard*(1-is_quasi_steady)+is_quasi_steady;
    % Number of cycles to run
    params.n_cycles = params.n_discard+1;
    % Time of one cycle
    params.t_cycle = t_cycle;
    % Total test time
    params.tf = params.t_cycle*params.n_cycles;    
    % Time span array
    params.tspan = [0 params.tf];                           
end

%% Set unavailable experimental/model data to NaN
data_fields_list = {"alpha_exp_cl","alpha_exp_cm","alpha_exp_cd","alpha_exp_cn","alpha_exp_cc","alpha_exp_ch",...
                    "delta_exp_cl","delta_exp_cm","delta_exp_cd","delta_exp_cn","delta_exp_cc","delta_exp_ch",...
                    "cl_exp","cm_exp","cd_exp","cn_exp","cc_exp","ch_exp",...
                    "time_exp_cl","time_exp_cm","time_exp_cd","time_exp_cn","time_exp_cc",...
                    "clt_exp","cmt_exp","cdt_exp","cnt_exp","cct_exp",...
                    "alpha_mod_cl","alpha_mod_cm","alpha_mod_cd","alpha_mod_cn","alpha_mod_cc","alpha_mod_ch",...
                    "delta_mod_cl","delta_mod_cm","delta_mod_cd","delta_mod_cn","delta_mod_cc","delta_mod_ch",...
                    "cl_mod","cm_mod","cd_mod","cn_mod","cc_mod","ch_mod",...
                    "time_mod_cl","time_mod_cm","time_mod_cd","time_mod_cn","time_mod_cc",...
                    "clt_mod","cmt_mod","cdt_mod","cnt_mod","cct_mod",...
                    "time_cycles","alpha_cycles","cl_cycles","cm_cycles","cd_cycles","cn_cycles","cc_cycles"};
for i=1:length(data_fields_list)
    field = data_fields_list{i};
    if ~isfield(data,field)
        data.(field) = NaN;
    end
end

%% Airfoil parameters
switch model
    case {"BL","BLT","BLTflap","BLgust","BLTgust","BLtvfs","BLTtvfs","BLhargen","BLThargen"}
        params = BL_airfoil_parameters(params,airfoil,M,U,b);
    case {"BLS","BLSLR"}
        params = BLS_airfoil_parameters(params,airfoil,M,U,b);
    case {"BLO","BLOgust"}
        params = BLO_airfoil_parameters(params,airfoil,M,U,b);  
    case "BLG"
        params = BLG_airfoil_parameters(params,airfoil,M,U,b);  
    otherwise
        error("Unknown model: " + model);
end

%% Indicial parameters
params = indicial_params(params,model,gust_ind);

%% State space matrices
params = SS_matrices(params,model,gust_ind);

%% Modifications for specific cases
params = modify_params(params,data,source,model,gust_ind,case_now);

end