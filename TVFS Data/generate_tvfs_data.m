clc
clear

% Experimental setup
tvfs_case = "05";
authors = {"CFD - Jose (2006)","Model - Jose (2006)"};
airfoil = "NACA0012";
M = 0.5;
lambda = 0.8;
k = 0.0;
k_U = 0.2;
a_0 = 1*pi/180;
a_1 = 0*pi/180;
psi_a = 0;
b = 0.3; % This is unknown, but the compressible flow inertial loads depend on it! 

a_inf = 340;
U = M*a_inf;
beta = sqrt(1-M^2);
U_0 = U;
U_1 = lambda*U;

% Available coefficients data
coefs = ["cl","cm"];

% Load data points
for i=1:length(coefs)
    % Load
    filename_exp = "exp_" + coefs(i) + "_" + tvfs_case + ".mat";
    filename_mod = "mod_" + coefs(i) + "_" + tvfs_case + ".mat";
    load(filename_exp); load(filename_mod); 
end

% Set data arrays
time_exp_cl = exp_cl_data(:,1);
time_mod_cl = mod_cl_data(:,1);
clt_exp = exp_cl_data(:,2);
clt_mod = mod_cl_data(:,2);
time_exp_cm = exp_cm_data(:,1);
time_mod_cm = mod_cm_data(:,1);
cmt_exp = exp_cm_data(:,2);
cmt_mod = mod_cm_data(:,2);

% Save file
new_filename = "tvfs_case_" + tvfs_case + ".mat";
save(new_filename,'authors','M','U','b','a_inf','beta','a_0','a_1','k','U_0','U_1','k_U','psi_a','airfoil','time_exp_cl','clt_exp','time_mod_cl','clt_mod','time_exp_cm','cmt_exp','time_mod_cm','cmt_mod');
