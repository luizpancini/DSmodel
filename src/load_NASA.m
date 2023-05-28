function [params,data] = load_NASA(frame)

% Load experiment condition variables
filepath = sprintf('../NASA Data/frame_%d.mat',frame);        
load(filepath,'authors','airfoil','a_inf','b','ah','U','M','beta','a_0','a_1','k','alpha_exp_cl','cl_exp','alpha_exp_cm','cm_exp','alpha_exp_cd','cd_exp','time_exp_cl','clt_exp','time_exp_cm','cmt_exp','time_exp_cd','cdt_exp');

% Set experimental data on struct
data = variables2struct(struct(),authors,alpha_exp_cl,alpha_exp_cm,alpha_exp_cd,cl_exp,cm_exp,cd_exp,time_exp_cl,time_exp_cm,time_exp_cd,clt_exp,cmt_exp,cdt_exp);

% Set all flow and test condition variables into the params struct
a_1h = 0; k_h = 0;
params = variables2struct(struct(),M,U,b,a_inf,beta,a_0,a_1,k,airfoil,ah,a_1h,k_h);

end

