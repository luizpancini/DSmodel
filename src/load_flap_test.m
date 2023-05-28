function [params,data] = load_flap_test(flap_case)

% Load experimental/modeling condition variables
filepath = sprintf('../Flap Data/flap_case_%02d.mat', flap_case);        
load(filepath,'authors','airfoil','b','a_inf','M','U','beta','ah','dh','k_f','d_0','d_1','a_0','a_1','k','delta_exp_cn','delta_exp_cm','delta_exp_ch','cn_exp','cm_exp','ch_exp','delta_mod_cn','delta_mod_cm','delta_mod_ch','cn_mod','cm_mod','ch_mod');

% Set experimental/modeling data on struct
data = variables2struct(struct(),authors,delta_exp_cn,delta_exp_cm,delta_exp_ch,cn_exp,cm_exp,ch_exp,delta_mod_cn,delta_mod_cm,delta_mod_ch,cn_mod,cm_mod,ch_mod);

% Set all flow and test condition variables into the params struct
a_1h = 0; k_h = 0;
params = variables2struct(struct(),M,U,b,ah,dh,a_inf,beta,a_0,a_1,k,d_0,d_1,k_f,airfoil,a_1h,k_h);

end

