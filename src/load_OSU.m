function [params,data] = load_OSU(OSU)

% Load experiment condition variables
filepath = sprintf('../OSU Data/mat_files/OSU_%d.mat',OSU);        
load(filepath,'authors','airfoil','b','a_inf','Re','M','U','beta','k','a_0','a_1','ah','time','alpha','cl','cd','cm','cn','cc','time_cycles','alpha_cycles','cl_cycles','cd_cycles','cm_cycles','cn_cycles','cc_cycles');

% Set experimental data on struct
data.authors = authors;
data.alpha_exp_cl = alpha;
data.alpha_exp_cm = alpha;
data.alpha_exp_cd = alpha;
data.alpha_exp_cn = alpha;
data.alpha_exp_cc = alpha;
data.cl_exp = cl;
data.cm_exp = cm;
data.cd_exp = cd;
data.cn_exp = cn;
data.cc_exp = cc;
data.time_exp_cl = time;
data.time_exp_cm = time;
data.time_exp_cd = time;
data.time_exp_cn = time;
data.time_exp_cc = time;
data.clt_exp = cl;
data.cmt_exp = cm;
data.cdt_exp = cd;
data.cnt_exp = cn;
data.cct_exp = cc;
data.time_cycles = time_cycles;
data.alpha_cycles = alpha_cycles;
data.cl_cycles = cl_cycles;
data.cm_cycles = cm_cycles;
data.cd_cycles = cd_cycles;
data.cn_cycles = cn_cycles;
data.cc_cycles = cc_cycles; 

% Set all flow and test condition variables into the params struct
a_1h = 0; k_h = 0;
params = variables2struct(struct(),M,U,b,a_inf,beta,a_0,a_1,k,airfoil,ah,a_1h,k_h);

end

