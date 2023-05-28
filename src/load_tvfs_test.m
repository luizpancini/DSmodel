function [params,data] = load_tvfs_test(tvfs_case)

% Load experimental/modeling condition variables
filepath = sprintf('../TVFS Data/tvfs_case_%02d.mat', tvfs_case);
warning('off')
load(filepath,'authors','M','U','b','a_inf','beta','a_0','a_1','k','U_0','U_1','k_U','psi_a','airfoil','lambda_U','ah','time_exp_cl','clt_exp','time_mod_cl','clt_mod','time_exp_cm','cmt_exp','time_mod_cm','cmt_mod');
warning('on')

% Set experimental/modeling data on struct
data = variables2struct(struct(),authors,time_exp_cl,clt_exp,time_mod_cl,clt_mod,time_exp_cm,cmt_exp,time_mod_cm,cmt_mod);

% Set all flow and test condition variables into the params struct
a_1h = 0; k_h = 0;
params = variables2struct(struct(),M,U,b,a_inf,beta,a_0,a_1,k,U_0,U_1,k_U,psi_a,airfoil,lambda_U,ah,a_1h,k_h);

end