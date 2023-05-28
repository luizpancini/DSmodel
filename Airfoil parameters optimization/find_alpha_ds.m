clc
clearvars -except reset_errors E_cn_vec_best E_cm_vec_best E_cc_vec_best
addpath('../functions')

%% Inputs
% Input options:
frame = 8306;                                          
GUD = 11012812;            
OSU = [372 373 374 384 385 386 396 397 398 408 409 410 422 432 434];           
% Select input option as "NASA", "GU", "OSU"
INPUTS.source = "OSU";
% Select model (BL or BLS, for Sheng et al.'s model)
INPUTS.model = "BL";
% For BLS model, choose linearized reattachment or not
INPUTS.lin_reat = 0;
% Plot option
plotting = 1;

%% Loop over cases
% Pre-process
switch INPUTS.source
    case "NASA"
        INPUTS.cases = frame;
        input_name = "Frame ";
    case "GU"
        INPUTS.cases = GUD;
        input_name = "GUD ";
    case "OSU"
        INPUTS.cases = OSU;
        input_name = "OSU run ";
end
if ~isrow(INPUTS.cases), INPUTS.cases = INPUTS.cases'; end
% if plotting, close all; end

% Initialize error sums for current run
E_cn = 0;
E_cm = 0;
E_cc = 0;

% Initialize error vectors for each case
N = length(INPUTS.cases);
E_cn_vec = zeros(N,1);
E_cm_vec = zeros(N,1);
E_cc_vec = zeros(N,1);

% Initialize best error vectors if non-existent
if ~exist('E_cn_vec_best','var')
    E_cn_vec_best = nan(N,1);
    E_cm_vec_best = nan(N,1);
    E_cc_vec_best = nan(N,1);
end

% Initialize containers
alpha_ds_exp_vec = nan(N,1);
alpha_ds_th_vec = nan(N,1);
alpha_ds_cc_vec = nan(N,1);
alpha_lag_ds_th_vec = nan(N,1);
alpha_lag_ds_cc_vec = nan(N,1);
r_ds_exp_vec = nan(N,1);
r_ds_th_vec = nan(N,1);
r_ds_cc_vec = nan(N,1);
RD_ds_th_vec = nan(N,1);
RD_ds_cc_vec = nan(N,1);
cn_re_exp_vec = nan(N,1);
alpha_re_exp_vec = nan(N,1);
cn_re_vec = nan(N,1);
alpha_re_mincn_vec = nan(N,1);
M_vec = nan(N,1);
% Loop
j = 0;
for case_now = INPUTS.cases
    
    j = j+1;
    
    %% Initialization
    [data,params,x0,y0,ODE_options] = initialize_BL_model(case_now,INPUTS);
    
    %% ODEs solver
    tic;
    switch INPUTS.model
        case "BL"
            [t,x,y,xdot] = BL_RKF45(params.tspan,x0,y0,params,ODE_options);
    end
    run_time = toc;
    str = input_name + num2str(case_now) + ", runtime: " + num2str(run_time) + " s"; fprintf(str);

    %% Get output variables
    switch INPUTS.model
        case "BL"
            outputs = BL_output_vars(x,y,t);
    end
    
    %% Interpolate data and find normalized errors
    interp = interp_data(INPUTS,data,params,outputs,1,E_cn_vec_best(j),E_cm_vec_best(j),E_cc_vec_best(j));
    
    %% Set errors for current case and increment error sums
    E_cn_vec(j) = interp.NRMSE.cn;
    E_cm_vec(j) = interp.NRMSE.cm;
    E_cc_vec(j) = interp.NRMSE.cc;
    E_cn = E_cn + E_cn_vec(j);
    E_cm = E_cm + E_cm_vec(j);
    E_cc = E_cc + E_cc_vec(j);
    
    %% Plots
    if plotting
        call_plotter(INPUTS,case_now,data,params,outputs,interp);
    end
    
    %% Model's dynamic stall data
    % Unpack case parameters and outputs
    t_cycle = params.t_cycle;
    a_0 = 180/pi*params.a_0;
    a_1 = 180/pi*params.a_1;
    b = params.b;
    U = params.U;
    k = params.k;
    M = params.M;
    r0 = params.r0;
    t = outputs.t;
    alpha = 180/pi*outputs.alpha;
    theta = outputs.theta;
    qR = outputs.qR;
    c_n = outputs.c_n;
    c_m = outputs.c_m;
    c_c = outputs.c_c;
    t_interp = interp.time;
    alpha_interp = interp.alpha;
    cc_interp = interp.cc;
    cn_interp = interp.cn;
    cm_interp = interp.cm;
    N_interp = length(t_interp);
    half_N_interp = round(N_interp/2);
    
    % Indices
    itt = find(t>t(end)-5/4*t_cycle,1,'first');     % Index of cycle angle -90° for last cycle
    ittf = find(t>t(end)-1/4*t_cycle,1,'first');    % Index of cycle angle 270° for last cycle
    lc_range = itt:ittf;                            % Last cycle's full index range
    lc_down_range = itt+round((ittf-itt)/2):ittf;   % Last cycle's downstroke index range      
    interp_down_range = half_N_interp:N_interp;
    
    % Vectors for last cycle
    alpha_lc = alpha(lc_range);
    cn_lc = c_n(lc_range);
    cm_lc = c_m(lc_range);
    cc_lc = c_c(lc_range);
    
    % Stall onset - theta = 1 criterion
%     i_ds = find(theta(lc_range)>1,1,'first')+itt-1;     % Index at which theta > 1
%     if isempty(i_ds), disp("No-stall case, disconsidered..."); continue; end
%     alpha_ds_th = alpha(i_ds);                          % Angle of incidence when theta > 1
%     r_ds_th = qR(i_ds)*r0;                              % Reduced pitch rate when theta > 1
%     alpha_lag_ds_th = 180/pi*x(9,i_ds);                 % Lagged angle of attack when theta > 1
%     RD_ds_th = x(14,i_ds);                              % Lagged capped reduced pitch rate ratio when theta > 1
%     
%     % Stall onset - maximum c_c criterion
%     [cc_max_mod,ind_max_cc_mod] = max(cc_lc);                       % Value and index of maximum c_c over the last cycle's range
%     ind_max_cc_mod = ind_max_cc_mod+itt-1;                          % Index of maximum c_c over the overall time range
%     time_max_cc_mod = (t(ind_max_cc_mod)-t(itt))/t_cycle*360-90;    % Time (cycle angle) of maximum c_c
%     alpha_ds_cc = alpha(ind_max_cc_mod);                            % Angle of incidence of maximum c_c
%     r_ds_cc = qR(ind_max_cc_mod)*r0;                                % Reduced pitch rate of maximum c_c 
%     alpha_lag_ds_cc = 180/pi*x(9,ind_max_cc_mod);                   % Lagged angle of attack of maximum c_c 
%     RD_ds_cc = x(14,ind_max_cc_mod);                                % Lagged capped reduced pitch rate ratio of maximum c_c 
    
%     % Reattachment angle - maximum c_m criterion
%     [~,ind_re_maxcm] = max(c_m(lc_down_range));                 % Index for reattachment angle according to max c_m criterion
%     ind_re_maxcm = ind_re_maxcm+itt-1+round((ittf-itt)/2);      % Index of max c_m over the overall time range
%     t_re_maxcm = (t(ind_re_maxcm)-t(itt))/t_cycle*360-90;       % Cycle angle when reattachment angle is reached (for max c_m criterion)
%     alpha_re_maxcm = alpha(ind_re_maxcm);                       % Angle of incidence when reattachment angle is reached (for max c_m criterion)
%     r_re_maxcm = qR(ind_re_maxcm)*r0;                           % Reduced pitch rate when reattachment angle is reached (for max c_m criterion)
%    
    % Experimental reattachment angle - minimum c_n criterion
%     ind_lmin_cn = interp_down_range(1) + find(islocalmin(cn_interp(interp_down_range))==1);             % Indices where c_n is a minimum in the downstroke
% %     figure;plot(t_interp,cn_interp,t_interp(ind_lmin_cn),cn_interp(ind_lmin_cn),'ko');grid
%     if isempty(ind_lmin_cn)
%         ind_re_exp = round(interp_down_range(1));                                % If there is no local minimum, then assume value at begin of downstroke
%     else
%         [cn_min_down,ind_cn_min_down] = min(cn_interp(interp_down_range));  % Absolute minimum c_n in the downstroke
%         ind_cn_min_down = ind_cn_min_down + interp_down_range(1);           % Index for absolute minimum c_n in the downstroke
%         ind_re_exp = ind_lmin_cn(ind_cn_min_down == ind_lmin_cn);           % If the absolute minimum is a local minimum, then it is the reattachment value
%         if isempty(ind_re_exp)                                              % But if not
%             % Then set reattachment value as the smallest of the local minima
%             [~,ind] = min(cn_interp(ind_lmin_cn));
%             ind_re_exp = ind_lmin_cn(ind);
%         end
%     end
%     cn_re_exp = cn_interp(ind_re_exp);                                      % Local minimum c_n in the reattachment
%     alpha_re_exp = alpha_interp(ind_re_exp);                                 % Angle of incidence when reattachment angle is reached (for min c_n criterion)
% %     hold on; plot(t_interp(ind_re_exp),cn_re_exp,'rd'); hold off
%     % Reattachment angle - minimum c_n criterion
%     ind_re_mincn = find(islocalmin(c_n(lc_down_range))==1);             % Index for reattachment angle according to min c_n criterion
%     if isempty(ind_re_mincn)
%         ind_re_mincn = lc_down_range(1);                                % If there is no local minimum, then assume value at begin of downstroke
%     else
%         ind_re_mincn = ind_re_mincn(end)+itt-1+round((ittf-itt)/2);     % Index of min c_n over the overall time range     
%     end    
%     cn_re = c_n(ind_re_mincn);                                 % Local minimum c_n in the reattachment
%     t_re_mincn = (t(ind_re_mincn)-t(itt))/t_cycle*360-90;       % Cycle angle when reattachment angle is reached (for min c_n criterion)
%     alpha_re_mincn = alpha(ind_re_mincn);                       % Angle of incidence when reattachment angle is reached (for min c_n criterion)
%     r_re_mincn = qR(ind_re_mincn)*r0;                           % Reduced pitch rate when reattachment angle is reached (for min c_n criterion)   
    
%     % Experimental extrema loads
%     [cn_max_exp,ind_max_cn_exp] = max(cn_interp);           % Value and index of maximum c_n over the last cycle's range
%     time_max_cn_exp = t_interp(ind_max_cn_exp);             % Time (cycle angle) of maximum c_n
%     alpha_cn_max_exp = alpha_interp(ind_max_cn_exp);     % Angle of incidence of maximum c_n
%     [cm_max_exp,ind_max_cm_exp] = min(cm_interp);           % Value and index of maximum (negative) c_m over the last cycle's range
%     time_max_cm_exp = t_interp(ind_max_cm_exp);             % Time (cycle angle) of maximum c_m
%     alpha_cm_max_exp = alpha_interp(ind_max_cm_exp);     % Angle of incidence of maximum c_m
%     % Model extrema loads
%     [cn_max_mod,ind] = max(cn_lc); alpha_cn_max_mod = alpha_lc(ind); % Value and angle of incidence of maximum c_n
%     [cm_max_mod,ind] = min(cm_lc); alpha_cm_max_mod = alpha_lc(ind); % Value and angle of incidence of maximum (negative) c_m
%     figure; plot(alpha_interp,cn_interp,'k-',alpha_cn_max_exp,cn_max_exp,'ko',alpha_lc,cn_lc,'b-',alpha_cn_max_mod,cn_max_mod,'bo'); grid
%     figure; plot(alpha_interp,cm_interp,'k-',alpha_cm_max_exp,cm_max_exp,'ko',alpha_lc,cm_lc,'b-',alpha_cm_max_mod,cm_max_mod,'bo'); grid
%     % Loads extrema
%     [cn_max,ind] = max(cn_lc); angle_cn_max = (t(ind+itt-1)-t(itt))/t_cycle*360-90; r_cn_max = qR(ind+itt-1)*r0;
%     [cm_max,ind] = max(-cm_lc); angle_cm_max = (t(ind+itt-1)-t(itt))/t_cycle*360-90; r_cm_max = qR(ind+itt-1)*r0;
%     [cc_max,ind] = max(cc_lc); angle_cc_max = (t(ind+itt-1)-t(itt))/t_cycle*360-90; r_cc_max = qR(ind+itt-1)*r0;
    
%     %% Experimental dynamic stall - per maximum c_c criterion
%     [cc_max_exp,ind_max_cc_exp] = max(interp.cc);                    % Index of maximum c_c over the last cycle's range
%     time_max_cc_exp = t_interp(ind_max_cc_exp);             % Time (cycle angle) of maximum c_c  
%     alpha_ds_exp = a_0 + a_1*sind(time_max_cc_exp);         % Angle of incidence of maximum c_c
%     alpha_dot_ds_exp = a_1*k*U/b*cosd(time_max_cc_exp);     % Pitch rate of maximum c_c
%     q_ds_exp = pi/180*2*alpha_dot_ds_exp*b/U;               % Nondimensional pitch rate of maximum c_c 
%     r_ds_exp = abs(q_ds_exp)/2;                             % Reduced pitch rate of maximum c_c 
    
%     figure; plot(alpha_interp,cc_interp,'k-',alpha_ds_exp,cc_max_exp,'ko',alpha(lc_range),cc_lc,'b-',alpha_ds_cc,cc_max_mod,'bo'); grid
    
    %% Set resulting values in container vectors
%     alpha_ds_exp_vec(j) = alpha_ds_exp;
%     alpha_ds_th_vec(j) = alpha_ds_th;
%     alpha_ds_cc_vec(j) = alpha_ds_cc;
%     alpha_lag_ds_th_vec(j) = alpha_lag_ds_th;
%     alpha_lag_ds_cc_vec(j) = alpha_lag_ds_cc;
%     r_ds_exp_vec(j) = r_ds_exp;
%     r_ds_th_vec(j) = r_ds_th;
%     r_ds_cc_vec(j) = r_ds_cc;
%     RD_ds_th_vec(j) = RD_ds_th;
%     RD_ds_cc_vec(j) = RD_ds_cc;
%     cn_re_exp_vec(j) = cn_re_exp;                                     
%     alpha_re_exp_vec(j) = alpha_re_exp;
%     cn_re_vec(j) = cn_re;
%     alpha_re_mincn_vec(j) = alpha_re_mincn;
    
    % Increment for mean Mach number calculation
    M_vec(j) = M;
end

%% Display total error sums and update best
if isnan(E_cn_vec_best(1))
    str = "Mean errors: cn = " + num2str(mean(E_cn_vec_best),'%.4f') + ", cm = " + num2str(mean(E_cm_vec_best),'%.4f') + ", cc = " + num2str(mean(E_cc_vec_best),'%.4f');
else
    str = "Mean errors: cn = " + num2str(mean(E_cn_vec_best),'%.4f') + "(" + num2str(100*(E_cn/sum(E_cn_vec_best)-1),'%.1f') + "%), cm = " + num2str(mean(E_cm_vec_best),'%.4f') + "(" + num2str(100*(E_cm/sum(E_cm_vec_best)-1),'%.1f') + "%), cc = " + num2str(mean(E_cc_vec_best),'%.4f') + "(" + num2str(100*(E_cn/sum(E_cn_vec_best)-1),'%.1f') + "%)";
end
disp(str);
if E_cn < sum(E_cn_vec_best) || isnan(E_cn_vec_best(1))
    E_cn_vec_best = E_cn_vec;
end
if E_cm < sum(E_cm_vec_best) || isnan(E_cm_vec_best(1))
    E_cm_vec_best = E_cm_vec;
end
if E_cc < sum(E_cc_vec_best) || isnan(E_cc_vec_best(1))
    E_cc_vec_best = E_cc_vec;
end

%% Get mean parameters
M_mean = mean(M_vec,'omitnan');
switch INPUTS.model
    case "BL"
        p = airfoil_parameters_BL(params,params.airfoil,M_mean,M_mean*params.a_inf,b);
end
alpha_ds0 = p.alpha_ds0; alpha_ss = p.alpha_ss; r0 = p.r0; Ta = p.Ta;

%% Plots
%%%%%%
r_f = 0.02;
x_dashed = [0 r0];
y_dashed = [14*pi/180 alpha_ds0+Ta*r0];
yy_dashed = 180/pi*[alpha_ss alpha_ds0+Ta*r0];
x2_dashed = [r0 r_f];
y2_dashed = 180/pi*[alpha_ds0+Ta*r0 alpha_ds0+Ta*r0+Ta*(r_f-r0)];
%%%%%%
set(0,'DefaultTextInterpreter','tex')
set(0,'DefaultLegendInterpreter','tex')
axes_size = 20;
lw = 0.75;
ms = 5;
figure1 = figure('InvertHardcopy','off','Color',[1 1 1]); figure1.PaperSize = [6 4.5]; figure1.Name = "AoA vs. r";
axes1 = axes('Parent',figure1,'FontSize',axes_size,'FontName','times new roman');
hold(axes1,'on'); grid on;
xlabel('r','FontWeight','normal','FontSize',axes_size);
ylabel('\alpha_{ds} [deg]','FontWeight','normal','FontSize',axes_size);
p0 = plot(nan,nan,'bo--',nan,nan,'ko',nan,nan,'ro','MarkerSize',ms,'Parent',axes1);
lgd1 = legend(p0,'Mod. (\theta = 1)','Exp. (max c_{c})','Mod. (max c_{c})');
set(lgd1,'Location','best','FontSize',12,'Box','on','AutoUpdate','off');
ylim([14 24]);
p1 = plot(r_ds_exp_vec,alpha_ds_exp_vec,'ko',r_ds_cc_vec,alpha_ds_cc_vec,'ro',r_ds_th_vec,alpha_ds_th_vec,'bo','MarkerSize',ms,'Parent',axes1);
p2 = plot(x_dashed,yy_dashed,'b--',x2_dashed,y2_dashed,'b--','LineWidth',lw,'Parent',axes1);
for i=1:N
    p3 = plot([r_ds_exp_vec(i) r_ds_cc_vec(i)],[alpha_ds_exp_vec(i) alpha_ds_cc_vec(i)],'k:','MarkerSize',ms,'Parent',axes1);
end
