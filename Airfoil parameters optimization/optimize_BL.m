function [E_sum,ps,E] = optimize_BL(p,INPUTS,varargin)

% Handle inputs
plotting = 0;       % TF to plot coefficients 
plot_events = 0;    % TF to plot DS events for check
dt_max = 1e-3;      % Maximum time step
if ~isempty(varargin)
    plotting = varargin{1};
    if length(varargin) > 1
        plot_events = varargin{2};
    end
    if length(varargin) > 2
        dt_max = varargin{3};
    end
end

%% Loop over cases
% Pre-process
if ~isrow(INPUTS.cases), INPUTS.cases = INPUTS.cases'; end
switch INPUTS.source
    case "NASA"
        INPUTS.base_name = "Frame ";
    case "GU"
        INPUTS.base_name = "GU run ";
    case "OSU"
        INPUTS.base_name = "OSU run ";
    case "gust_test"
        INPUTS.base_name = "Gust test case ";
    case "flap_test"
        INPUTS.base_name = "Flap test case ";
    case "tvfs_test"
        INPUTS.base_name = "TVFS test case ";    
end

% Initialize error sums and vectors
n_cases = length(INPUTS.cases);
E_cn_sum = 0;           E_cn_vec = nan(n_cases,1);
E_cm_sum = 0;           E_cm_vec = nan(n_cases,1);
E_cc_sum = 0;           E_cc_vec = nan(n_cases,1);
E_alpha_ds_sum = 0;     E_alpha_ds_vec = nan(n_cases,1);
E_alpha_re_sum = 0;     E_alpha_re_vec = nan(n_cases,1);
E_alpha_cnmax_sum = 0;  E_alpha_cnmax_vec = nan(n_cases,1);
E_alpha_cmmax_sum = 0;  E_alpha_cmmax_vec = nan(n_cases,1);
E_cn_re_sum = 0;        E_cn_re_vec = nan(n_cases,1);
E_cnmax_sum = 0;        E_cnmax_vec = nan(n_cases,1);
E_cmmax_sum = 0;        E_cmmax_vec = nan(n_cases,1);
E_ccmax_sum = 0;        E_ccmax_vec = nan(n_cases,1);
E_sum = 0;

% Loop
j = 0;
for case_now = INPUTS.cases
    
    j = j+1;
    
    %% Initialization
    % Select number of cycles to discard
    if INPUTS.source == "OSU" && (408 <= case_now && case_now <= 419)
        n_discard = 2;
    else
        n_discard = 1;
    end
    [data,params,x0,y0,ODE_options] = initialize_model(case_now,INPUTS,n_discard,dt_max,0);
    % Override with input parameter values
    [params,ps] = override_BL_params(p,params,INPUTS.model);
    
    %% ODEs solver
    switch INPUTS.model
        case "BL"
            outputs = BL_RKF45(params.tspan,x0,y0,params,ODE_options);
    end
        
    %% Interpolate data and find normalized errors
    interp = interp_data(case_now,INPUTS,data,params,outputs,0);
    % Errors
    E_cn = interp.NRMSE.cn; E_cn_vec(j) = E_cn;
    E_cm = interp.NRMSE.cm; E_cm_vec(j) = E_cm;
    E_cc = interp.NRMSE.cc; E_cc_vec(j) = E_cc;
    
    %% Plots
    if plotting
        call_plotter(INPUTS,case_now,data,params,outputs,interp);
    end
    
    %% Errors on dynamic stall events
    % Initialize errors
    E_alpha_ds = 0;
    E_alpha_re = 0;
    E_cn_re = 0;
    E_cn_max = 0; E_alpha_cn_max = 0;
    E_cm_max = 0; E_alpha_cm_max = 0;
    E_cc_max = 0;
    
    if INPUTS.events_weights.ds > 0 || INPUTS.events_weights.re > 0 || INPUTS.events_weights.le > 0
        % Unpack case parameters and outputs
        t_cycle = params.t_cycle;
        a_1 = 180/pi*params.a_1;
        t = outputs.t;
        alpha = 180/pi*outputs.alpha;
        c_n = outputs.c_n;
        c_m = outputs.c_m;
        c_c = outputs.c_c;
        t_ref = interp.time_ref;
        alpha_ref = interp.alpha_ref;
        cc_ref = interp.cc_ref;
        cn_ref = interp.cn_ref;
        cm_ref = interp.cm_ref;
        N_interp = length(t_ref);
        half_N_interp = round(N_interp/2);
        % Indices
        itt = find(t>t(end)-5/4*t_cycle,1,'first');             % Index of cycle angle -90° for last cycle
        ittf = find(t>t(end)-1/4*t_cycle,1,'first');            % Index of cycle angle 270° for last cycle
        lc_range = itt:ittf;                                    % Last cycle's full index range
        lc_down_range = itt+round((ittf-itt)/2):ittf;           % Last cycle's downstroke index range
        interp_down_range = half_N_interp:N_interp;             % Index range for interpolated experimental values
        % Vectors for last cycle
        alpha_lc = alpha(lc_range);
        cn_lc = real(c_n(lc_range));
        cm_lc = real(c_m(lc_range));
        cc_lc = real(c_c(lc_range));
        
        %%%%%% Stall onset prediction %%%%%%
        if INPUTS.events_weights.ds > 0
            % Experimental stall onset - maximum c_c criterion
            [cc_max_exp,ind_max_cc_exp] = max(cc_ref);              % Value and index of maximum c_c over the last cycle's range
            alpha_ds_exp = alpha_ref(ind_max_cc_exp);               % Angle of incidence of maximum c_c
            % Model stall onset - maximum c_c criterion
            [cc_max_mod,ind_max_cc_mod] = max(cc_lc);               % Value and index of maximum c_c over the last cycle's range
            alpha_ds_mod = alpha_lc(ind_max_cc_mod);                % Angle of incidence of maximum c_c
            % Difference in angle of stall onset normalized by AoA range
            E_alpha_ds = abs(alpha_ds_exp-alpha_ds_mod)/a_1;
            % Plot for check
            if plot_events
                figure; plot(alpha_ref,cc_ref,'k-',alpha_ds_exp,cc_max_exp,'ko',alpha_lc,cc_lc,'b-',alpha_ds_mod,cc_max_mod,'bo');grid;xlabel('AoA [deg]');ylabel('c_c');
            end
        end
        
        %%%%%% Reattachment prediction %%%%%%
        if INPUTS.events_weights.re > 0
            % Experimental reattachment angle - minimum c_n criterion
            ind_lmin_cn = interp_down_range(1) + find(islocalmin(cn_ref(interp_down_range))==1);    % Indices where c_n is a minimum in the downstroke
            if isempty(ind_lmin_cn)
                ind_re_exp = round(interp_down_range(1));                                           % If there is no local minimum, then assume value at begin of downstroke
            else
                [~,ind_cn_min_down] = min(cn_ref(interp_down_range));       % Absolute minimum c_n in the downstroke
                ind_cn_min_down = ind_cn_min_down + interp_down_range(1);   % Index for absolute minimum c_n in the downstroke
                ind_re_exp = ind_lmin_cn(ind_cn_min_down == ind_lmin_cn);   % If the absolute minimum is a local minimum, then it is the reattachment value
                if isempty(ind_re_exp)                                      % But if not
                    % Then set reattachment value as the smallest of the local minima
                    [~,ind] = min(cn_ref(ind_lmin_cn));
                    ind_re_exp = ind_lmin_cn(ind);
                end
                % But if the absolute minimum close to the end of the
                % downstroke, then it is probably a numerical interpolation
                % oscillation, so check other available options
                alpha_re_exp = alpha_ref(ind_re_exp);
                while abs(alpha_ref(end)-alpha_re_exp) < rad2deg(0.2*params.a_1) && length(ind_lmin_cn) > 1
                    % Remove that absolute minimum
                    [~,ind] = min(cn_ref(ind_lmin_cn));
                    ind_lmin_cn(ind) = []; 
                    % Set as "second absolute minimum"
                    [~,ind] = min(cn_ref(ind_lmin_cn));
                    ind_re_exp = ind_lmin_cn(ind);
                    alpha_re_exp = alpha_ref(ind_re_exp);
                end
            end         
            % If the angle of reattachment is identified as very close to
            % the end of the downstroke, then it is probably a false one,
            % set as the value at begin of downstroke
            if abs(alpha_ref(end)-alpha_ref(ind_re_exp)) < 0.5
                ind_re_exp = round(interp_down_range(1));
            end
            alpha_re_exp = alpha_ref(ind_re_exp);                               % Angle of incidence when reattachment angle is reached (for min c_n criterion)
            cn_re_exp = cn_ref(ind_re_exp);                                     % Local minimum c_n in the reattachment
            % Model reattachment angle - minimum c_n criterion
            ind_re_mod = find(islocalmin(real(c_n(lc_down_range))));                  % Index for reattachment angle according to min c_n criterion
            if isempty(ind_re_mod)
                ind_re_mod = lc_down_range(1);                                  % If there is no local minimum, then assume value at begin of downstroke
            else
                ind_re_mod = ind_re_mod(end)+itt-1+round((ittf-itt)/2);         % Index of min c_n over the overall time range
            end
            cn_re_mod = c_n(ind_re_mod);                                        % Local minimum c_n in the reattachment
            alpha_re_mod = alpha(ind_re_mod);                                   % Angle of incidence when reattachment angle is reached (for min c_n criterion)
            % Difference in angle of reattachment normalized by AoA range
            E_alpha_re = abs(alpha_re_exp-alpha_re_mod)/a_1;
            % Difference in cn at reattachment normalized by cn span
            E_cn_re = abs(cn_re_exp-cn_re_mod)/(max(cn_ref)-min(cn_ref));
            % Plot for check
            if plot_events
                figure;plot(alpha_ref,cn_ref,'k-',alpha_re_exp,cn_re_exp,'ko',alpha_lc,cn_lc,'b-',alpha_re_mod,cn_re_mod,'bo');grid;xlabel('AoA [deg]');ylabel('c_n'); 
            end
        end
        
        %%%%%% Loads extrema prediction %%%%%% 
        if INPUTS.events_weights.le > 0
            % Experimental extrema loads
            [cn_max_exp,ind_max_cn_exp] = max(cn_ref);                          % Value and index of maximum c_n
            alpha_cn_max_exp =alpha_ref(ind_max_cn_exp);                        % Angle of incidence of maximum c_n
            [cm_max_exp,ind_max_cm_exp] = min(cm_ref);                          % Value and index of maximum (negative) c_m
            alpha_cm_max_exp = alpha_ref(ind_max_cm_exp);                       % Angle of incidence of maximum c_m
            % Model extrema loads
            [cn_max_mod,ind] = max(cn_lc); alpha_cn_max_mod = alpha_lc(ind);    % Value and angle of incidence of maximum c_n
            [cm_max_mod,ind] = min(cm_lc); alpha_cm_max_mod = alpha_lc(ind);    % Value and angle of incidence of maximum (negative) c_m
            % Difference in angles of loads extrema normalized by AoA range
            E_alpha_cn_max = abs(alpha_cn_max_exp-alpha_cn_max_mod)/a_1;
            E_alpha_cm_max = abs(alpha_cm_max_exp-alpha_cm_max_mod)/a_1;
            % Difference in loads extrema normalized by loads range
            E_cn_max = abs(cn_max_exp-cn_max_mod)/(cn_max_exp-min(cn_ref));
            E_cm_max = abs((cm_max_exp-cm_max_mod)/(cm_max_exp-max(cm_ref)));
            E_cc_max = abs((cc_max_exp-cc_max_mod)/(cc_max_exp-min(cc_ref)));
            % Plot for check
            if plot_events
                figure; plot(alpha_ref,cn_ref,'k-',alpha_cn_max_exp,cn_max_exp,'ko',alpha_lc,cn_lc,'b-',alpha_cn_max_mod,cn_max_mod,'bo');grid;xlabel('AoA [deg]');ylabel('c_n');
                figure; plot(alpha_ref,cm_ref,'k-',alpha_cm_max_exp,cm_max_exp,'ko',alpha_lc,cm_lc,'b-',alpha_cm_max_mod,cm_max_mod,'bo');grid;xlabel('AoA [deg]');ylabel('c_m');
            end
        end
    end
    % Set vectors
    E_alpha_ds_vec(j) = E_alpha_ds;
    E_alpha_re_vec(j) = E_alpha_re;
    E_cn_re_vec(j) = E_cn_re;
    E_cnmax_vec(j) = E_cn_max; 
    E_alpha_cnmax_vec(j) = E_alpha_cn_max;
    E_cmmax_vec(j) = E_cm_max;
    E_alpha_cmmax_vec(j) = E_alpha_cm_max;
    E_ccmax_vec(j) = E_cc_max;
      
    %% Increment error sums
    E_cn_sum = E_cn_sum + E_cn;
    E_cm_sum = E_cm_sum + E_cm;
    E_cc_sum = E_cc_sum + E_cc;
    E_alpha_ds_sum = E_alpha_ds_sum + E_alpha_ds;
    E_alpha_re_sum = E_alpha_re_sum + E_alpha_re;
    E_cn_re_sum = E_cn_re_sum + E_cn_re;
    E_alpha_cnmax_sum = E_alpha_cnmax_sum + E_alpha_cn_max;
    E_alpha_cmmax_sum = E_alpha_cmmax_sum + E_alpha_cm_max;
    E_cnmax_sum = E_cnmax_sum + E_cn_max;
    E_cmmax_sum = E_cmmax_sum + E_cm_max;
    E_ccmax_sum = E_ccmax_sum + E_cc_max;
    
    % Normalized cumulative error sum
    E_sum = E_sum + INPUTS.cases_weights(j)*(INPUTS.coefs_weights.cn*E_cn + INPUTS.coefs_weights.cm*E_cm + INPUTS.coefs_weights.cc*E_cc + ...
                                      INPUTS.events_weights.ds*E_alpha_ds + INPUTS.events_weights.re*E_alpha_re + INPUTS.events_weights.re*E_cn_re + ...
                                      INPUTS.events_weights.le*E_cn_max + INPUTS.events_weights.le*E_cm_max + INPUTS.events_weights.le*E_cc_max + ...
                                      INPUTS.events_weights.le*E_alpha_cn_max + INPUTS.events_weights.le*E_alpha_cm_max)/...
                                      (INPUTS.coefs_weights.cn + INPUTS.coefs_weights.cm + INPUTS.coefs_weights.cc + ...
                                      INPUTS.events_weights.ds + 2*INPUTS.events_weights.re + 5*INPUTS.events_weights.le);
end

% Pack error data into structure
E = variables2struct(struct,E_cn_sum,E_cm_sum,E_cc_sum,E_alpha_ds_sum,E_alpha_re_sum,E_alpha_cnmax_sum,E_alpha_cmmax_sum,E_cn_re_sum,E_cnmax_sum,E_cmmax_sum,E_ccmax_sum,E_sum,E_cn_vec,E_cm_vec,E_cc_vec,E_alpha_ds_vec,E_alpha_re_vec,E_alpha_cnmax_vec,E_alpha_cmmax_vec,E_cn_re_vec,E_cnmax_vec,E_cmmax_vec,E_ccmax_vec);
    
end