function [time,alpha,cn_interp,cm_interp,cc_interp,cl_interp,cd_interp] = interp_coefs_NASA(data,a_0,a_1,N_interp,interp_method)

% Unpack
time_exp_cl = data.time_exp_cl;
time_exp_cm = data.time_exp_cm;
time_exp_cd = data.time_exp_cd;
clt_exp = data.clt_exp;
cmt_exp = data.cmt_exp;
cdt_exp = data.cdt_exp;

%% Cycle angle (time) and pitch angle arrays
time = linspace(-90,270,N_interp);
alpha = a_0+a_1*sind(time);

%% Interpolate coefficients
% Interpolation method
cl_interp = interp1(time_exp_cl,clt_exp,time,interp_method);
cm_interp = interp1(time_exp_cm,cmt_exp,time,interp_method);
cd_interp = interp1(time_exp_cd,cdt_exp,time,interp_method);
cn_interp = cl_interp.*cos(alpha)+cd_interp.*sin(alpha);
cc_interp = cl_interp.*sin(alpha)-cd_interp.*cos(alpha);

end