function [y,alpha_re_i,cn_re_i,cm_re_i,cc_re_i,cn_re_f,cm_re_f,cc_re_f,tri,trf] = BLS_outputs(t,x,tv0,U,b,beta,k,a_0,a_1,M,K0,K1,K2,kappa,T_I,T_M,b1,b2,b3,A1,A2,A3,eta,E0,Df,c_m0,c_n_alpha,TvL,r0,alpha_ds0,alpha_ss,alpha1_0,S1,S2,Tf0,f0,fb,B1,B2,alpha_0L,c_d0,alpha_min0,Tr,r_dcn_da_reat,r_dcm_da_reat,r_dcc_da_reat,alpha_re_i,cn_re_i,cm_re_i,cc_re_i,cn_re_f,cm_re_f,cc_re_f,tri,trf,theta_max,lin_reat)

%% Kinematics
[alpha,~,~,~,q,~,~,~,~,~,R,~] = BLO_kinematics(t,U,b,k,a_0,a_1,r0);

%% Stall onset criterion
[so_state,so_lim,theta] = BLS_stall_onset(x(7),R,alpha_ss,alpha_ds0);        

%% Unsteady breakpoint angle
alpha1 = BLS_alpha_brk(alpha1_0,alpha_ss,alpha_ds0,R,alpha,q);

%% Time of stall onset
tau_v = max([0 t-tv0]);

%% Separation points
[f,f_prime] = BLS_sep_points(so_state,f0,fb,alpha,alpha1,alpha1_0,S1,S2);

%% Tf variation
Tf = BLS_time_constants(Tf0,theta,tau_v,TvL,alpha,q,x(8),fb);

%% Airload coefficients
% c_n
[c_n,c_nC,c_nf,c_nI,c_nv,alpha_E] = BLS_cn_coeff(x,U,b,M,beta,A1,A2,A3,b1,b2,b3,T_I,alpha,q,c_n_alpha,f,x(8),B1,alpha_0L,tau_v,TvL);
% c_m
[c_m,~,c_mf,c_mI,c_mv,dCP] = BLS_cm_coeff(x,c_nf,c_nv,M,alpha,q,T_I,K0,K1,K2,kappa,tau_v,TvL,c_m0,U,b,beta,T_M,B2,f_prime);
% c_c
c_c = BLS_cc_coeff(alpha,so_state,so_lim,x(8),theta,c_n_alpha,alpha_E,E0,Df,alpha_0L,c_d0,f_prime,eta);
% Override with linear approximation for reattachment phase
if lin_reat == 1    
    [c_n,c_m,c_c,alpha_re_i,cn_re_i,cm_re_i,cc_re_i,cn_re_f,cm_re_f,cc_re_f,tri,trf] = BLS_reattach(t,b,U,alpha,q,tau_v,TvL,a_0,a_1,c_n_alpha,alpha_min0,Tr,r_dcn_da_reat,r_dcm_da_reat,r_dcc_da_reat,c_n,c_m,c_c,alpha_re_i,cn_re_i,cm_re_i,cc_re_i,cn_re_f,cm_re_f,cc_re_f,tri,trf,theta_max);
end

%% Outputs
y = [alpha; c_n; c_m; c_c; f; f_prime; c_nf; c_nI; c_nv; c_mf; c_mv; c_mI; tau_v; dCP; so_lim; q; alpha1; R; Tf; c_nC];

end