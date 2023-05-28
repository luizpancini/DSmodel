function [time,alpha,cn_interp,cm_interp,cc_interp,cl_interp,cd_interp,interp,itt,ittf,range] = interp_coefs_GU(data,a_0,a_1,t,t_cycle,interp,N_interp,interp_method)

% Unpack
alpha_exp = data.alpha_exp_cl;
time_exp = data.time_exp_cl;
cn_exp = data.cn_exp;
cm_exp = data.cm_exp;
cc_exp = data.cc_exp;
cl_exp = data.cl_exp;
cd_exp = data.cd_exp;

%% Find time shift
alpha_init = alpha_exp(1)*pi/180;           % Pitch angle at which the experiment begins
t_init = asin((alpha_init-a_0)/a_1)*180/pi; % Cycle angle [deg] at which the data starts
delta_t = -t_init/360;                      % According time shift (fraction of t_cycle)

%% Find new indices and update interpolation data structure
itt = find(t>t(end)-(1+delta_t)*t_cycle,1,'first');
ittf = min([length(t),find(t>t(end)-(delta_t)*t_cycle,1,'first')]);
range = itt:ittf;
interp.itt = itt;
interp.ittf = ittf;
interp.range = itt:ittf;

%% Cycle angle (time) and AoA vectors
time = linspace(-90,270,N_interp);
alpha = a_0+a_1*sind(time);

%% Interpolate coefficients
cn_interp = interp1(time_exp,cn_exp,time,interp_method);
cm_interp = interp1(time_exp,cm_exp,time,interp_method);
cc_interp = interp1(time_exp,cc_exp,time,interp_method);
cl_interp = interp1(time_exp,cl_exp,time,interp_method);
cd_interp = interp1(time_exp,cd_exp,time,interp_method);

end