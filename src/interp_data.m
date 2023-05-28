function interp = interp_data(case_now,INPUTS,data,params,outputs,varargin)

%% Handle inputs
% Unpack parameters and outputs
a_0 = params.a_0;
a_1 = params.a_1;
t = outputs.t;
alpha = outputs.alpha;
c_n = outputs.c_n;
c_m = outputs.c_m;
c_c = outputs.c_c;
source = INPUTS.source;
if source == "other" && any(params.k_h)
    alpha = outputs.alpha_plunge;
    a_1 = params.a_1h;
end

% Check if plot style was input
if isfield(params,'time_plot_style')
    time_plot_style = params.time_plot_style;
else
    time_plot_style = "";
end

% Find time indices of oscillatory cases
if isfield(params,'t_cycle')
    t_cycle = params.t_cycle;
    n_discard = params.n_discard;
    if contains(source,"tvfs") || time_plot_style == "tvfs"
        delta_cycle = 0; 
    else
        delta_cycle = 1/4; 
    end
    iti = find(t>n_discard*t_cycle,1,'first');                                      % First index after discarded time
    itt = find(t>t(end)-(1+delta_cycle)*t_cycle,1,'first');                         % Index of begin of last cycle
    ittf = min([length(t),find(t>t(end)-delta_cycle*t_cycle,1,'first')]);           % Index of end of last cycle
    range = itt:ittf;                                                               % Range of indices in last cycle
else
    iti = nan; itt = nan; ittf = nan; range = nan;
end

% Create interpolation data structure
interp = variables2struct(struct(),iti,itt,ittf,range);

% Check if general harmonic test contains data from references
if contains(source,"hargen") && isfield(data,'source'), source = data.source; end

% Check data availability for interpolation
if ~ismember(source,{'NASA','GU','OSU','other'}) 
    time_ref = nan; alpha_ref = nan; cn_ref = nan; cm_ref = nan; cc_ref = nan; cl_ref = nan; cd_ref = nan;
    interp = variables2struct(interp,time_ref,alpha_ref,cn_ref,cm_ref,cc_ref,cl_ref,cd_ref);
    return;
end

% Check optional inputs
disp_info = 1; error_standard = 'NRMSE'; E_cn_best = nan; 
j = 0;
while ~isempty(varargin)
    j = j+1;
    switch j
        case 1
            disp_info = varargin{1};
            varargin(1) = [];
        case 2
            error_standard = varargin{1};
            varargin(1) = [];
        case 3
            E_cn_best = varargin{1};
            E_cm_best = varargin{2};
            E_cc_best = varargin{3};
            varargin(1:3) = [];            
    end
end

%% Interpolation of coefficients
% Interpolation method
N_interp = 1e4+1;            % Size of interpolated arrays - must be kept constant for comparison between any two sets, because the NRMSE depends on it
interp_method = 'pchip';     % Interpolation method (pchip is best for OSU data)

% Interpolate reference data
switch source
    case "NASA"
        [time_ref_interp,alpha_ref_interp,cn_ref_interp,cm_ref_interp,cc_ref_interp,cl_ref_interp,cd_ref_interp] = interp_coefs_NASA(data,a_0,a_1,N_interp,interp_method);
    case "GU"
        [time_ref_interp,alpha_ref_interp,cn_ref_interp,cm_ref_interp,cc_ref_interp,cl_ref_interp,cd_ref_interp,interp,itt,ittf,range] = interp_coefs_GU(data,a_0,a_1,t,t_cycle,interp,N_interp,interp_method);
    case "OSU"
        [time_ref_interp,alpha_ref_interp,cn_ref_interp,cm_ref_interp,cc_ref_interp,cl_ref_interp,cd_ref_interp] = interp_coefs_OSU(data,N_interp,interp_method);
    case "other"
        [time_ref_interp,alpha_ref_interp,cn_ref_interp,cm_ref_interp,cc_ref_interp,cl_ref_interp,cd_ref_interp,interp] = interp_coefs_other(case_now,data,a_0,a_1,t,t_cycle,interp,N_interp,interp_method);
end

% Interpolate model's results
time_mod_interp = linspace(t(itt),t(ittf),N_interp);                                % Linearly interpolated time array
alpha_mod_interp = interp1(t(range),alpha(range),time_mod_interp,interp_method);    % Interpolated pitch angle
cn_mod_interp = interp1(t(range),c_n(range),time_mod_interp,interp_method);         % Interpolated c_n
cm_mod_interp = interp1(t(range),c_m(range),time_mod_interp,interp_method);         % Interpolated c_m
cc_mod_interp = interp1(t(range),c_c(range),time_mod_interp,interp_method);         % Interpolated c_c
% cl_mod_interp = cn_mod_interp.*cos(alpha_mod_interp)+cc_mod_interp.*sin(alpha_mod_interp);
% cd_mod_interp = cn_mod_interp.*sin(alpha_mod_interp)-cc_mod_interp.*cos(alpha_mod_interp);

%% Error calculation
switch error_standard
    case 'NRMSE'
        % Normalized Root Mean Squared Error
        [cn_error,cm_error,cc_error] = NRMSE_calculator(cn_mod_interp,cm_mod_interp,cc_mod_interp,cn_ref_interp,cm_ref_interp,cc_ref_interp);
    case 'RMSE'
        % Root Mean Squared Error
        [cn_error,cm_error,cc_error] = RMSE_calculator(cn_mod_interp,cm_mod_interp,cc_mod_interp,cn_ref_interp,cm_ref_interp,cc_ref_interp);
end

%% Display information
if disp_info
    if isnan(E_cn_best)
        str = ", Errors: cn = " + num2str(cn_error,'%.4f') + ", cm = " + num2str(cm_error,'%.4f') + ", cc = " + num2str(cc_error,'%.4f');
    else
        str = ", Errors: cn = " + num2str(cn_error,'%.4f') + "(" + num2str(100*(cn_error/E_cn_best-1),'%.1f') + "%%), cm = " + num2str(cm_error,'%.4f') + "(" + num2str(100*(cm_error/E_cm_best-1),'%.1f') + "%%), cc = " + num2str(cc_error,'%.4f') + "(" + num2str(100*(cc_error/E_cc_best-1),'%.1f') + "%%)";
    end
    fprintf(str);
end

%% Save interpolation data to struct
interp.time_mod = time_mod_interp;
interp.alpha_mod = alpha_mod_interp;
interp.cn_mod = cn_mod_interp;
interp.cm_mod = cm_mod_interp;
interp.cc_mod = cc_mod_interp;
interp.time_ref = time_ref_interp;
interp.alpha_ref = alpha_ref_interp;
interp.cn_ref = cn_ref_interp;
interp.cm_ref = cm_ref_interp;
interp.cc_ref = cc_ref_interp;
interp.cl_ref = cl_ref_interp;
interp.cd_ref = cd_ref_interp;
interp.error.cn = cn_error;
interp.error.cm = cm_error;
interp.error.cc = cc_error;

end