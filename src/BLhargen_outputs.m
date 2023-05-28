function [y,t_ub,t_db,r,qR,theta,alpha_lag,dalpha1_n,fprime_n,fprime_m,fprime_c,fprime_n_db,fprime_m_db,fprime_c_db,RD_ub,RD_db,theta_min,theta_max,RD_m,qR_max,S,P,upstroke,T_flag_downstroke,c_n_alpha,tv0,f_diff_tv0,qR_tv0,R_tv0,RD_tv0,Tv_tv0,theta_tv0,upstroke_tv0,f_diff_tv0_2,RD_tv0_2,upstroke_tv0_2,U,beta] = ...
         BLhargen_outputs(t,t_i,x,dt,tv0,t_ub,t_db,r_i,qR_i,theta_i,alpha_lag_i,upstroke_i,S_i,P_i,T_flag_downstroke_i,dalpha1_n_i,dalpha1_db_i,fprime_n_i,fprime_m_i,fprime_c_i,fprime_n_db_i,fprime_m_db_i,fprime_c_db_i,RD_ub_i,RD_db_i,theta_max_i,theta_min_i,RD_m_i,qR_max_i,c_n_alpha_i,qR_tv0,R_tv0,RD_tv0,f_diff_tv0,theta_tv0,upstroke_tv0,Tv_tv0,f_diff_tv0_2,RD_tv0_2,upstroke_tv0_2,airfoil,b,ah,dh,a_inf,a_0,a_1,k,a_1h,k_h,psi_ha,U_0,U_1,k_U,psi_Ua,d_0,d_1,k_f,psi_fa,lambda_g0,AG,bG,G,Cg,b5,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,alpha_0L,alpha1_0n,alpha1_0m,alpha1_0c,alpha_ds0,alpha_ss,beta_Sig1n,beta_Sig1c,beta_Sig2n,beta_S2n_lpr,beta_S2c_lpr,beta_S1n_u,beta_S1m_u,beta_S1c_u,beta_S1n_d,beta_S1m_d,beta_S1c_d,beta_S2n_u,beta_S2m_u,beta_S2c_u,beta_S2n_d,beta_S2m_d,beta_S2c_d,gamma_LS,delta_alpha_0,delta_alpha_1,eps_fn,eps_fm,eps_fh,eta,kappa_0,kappa_1,kappa_2,kappa_3,lambda_1,lambda_2,mu_v2,nu_1,nu_2,nu_3,nu_4,nu_5,chi_u,chi_d,xi,zeta_a,c_d0,c_m0,c_n_alpha,d_cc,d_cm,E0,E1,f0_n,f0_m,f0_c,fb_n,fb_m,fb_c,g_v,g_v2,K0,K1,K2,r0,S1_n,S1_m,S1_c,S2_n,S2_m,S2_c,Ta,Tf,Tv,Tv2,Vn1,Vn2,Vn3,Vm,Vc,z_ccd,z_ccu,z_cm,wg_fun,wgdot_fun,wgddot_fun)

%% Kinematics
[U,Udot,M,Mdot,beta,w_QS,Uc,alpha_QS,alphadot_QS,alpha,alphadot,alphaddot,hdot,hddot,alpha_plunge,delta,deltadot,deltaddot,w_tqc_pp,wdot_tqc_pp,cnawdot_tqc,q,q_QS,r,rdot,qR,R,w_tqc_f,wdot_tqc_f,Ucm_QS,Ucm_QS_dot,c_n_alpha] = ...
BLhargen_kinematics(t,dt,airfoil,b,ah,a_inf,a_0,a_1,k,a_1h,k_h,psi_ha,U_0,U_1,k_U,psi_Ua,d_0,d_1,k_f,psi_fa,T10,T11,r0,r_i,c_n_alpha_i);
% Gust-induced kinematics
[wg,vg] = gust_kin_comp(t,U,U_0,lambda_g0,wg_fun);

%% States
alpha_lag = atan((x(1)+x(7))/(Uc+vg));
f2prime_n = x(2);
f2prime_m = x(3);
f2prime_c = x(4);
RD = x(5);
RD_theta = x(6);
gust_states = x(23:end);

%% Stall onset criterion, stall and motion qualifiers
[alpha_cr,theta,theta_max,theta_min,del_RD_acc,RD_m,qR_max,R_max,S,P,T,upstroke,downstroke_beginning,in_stall,T_flag_downstroke,t_ub,t_db,T_ub,T_db,RD_ub,RD_db,T_s] = BL_stall_onset(t,t_ub,t_db,tv0,RD_ub_i,RD_db_i,RD_tv0,theta_i,alpha_lag_i,qR_i,upstroke_i,S_i,P_i,T_flag_downstroke_i,theta_max_i,theta_min_i,RD_m_i,qR_max_i,alpha_lag,q_QS,qR,R,RD,RD_theta,rdot,alpha_ds0,alpha_ss,gamma_LS,Tf,zeta_a);

%% Unsteady breakpoint of separation angles
[alpha1_n,dalpha1_n,alpha1_m,dalpha1_m,alpha1_c,dalpha1_c] = BL_alpha_brk(qR,R,RD,RD_theta,del_RD_acc,RD_m,theta,theta_min,theta_max,R_max,S,P,T,T_s,T_ub,T_db,upstroke,downstroke_beginning,dalpha1_n_i,dalpha1_db_i,alpha1_0n,alpha1_0m,alpha1_0c,alpha_ds0,alpha_ss,delta_alpha_0,delta_alpha_1,d_cm,d_cc,z_cm,z_ccd,z_ccu);
  
%% Separation points
[f_n,f_m,f_c,fprime_n,fprime_m,fprime_c,~,fprime_n_db,fprime_m_db,fprime_c_db] = BL_sep_points(alpha_QS,alpha_lag,alpha1_n,alpha1_m,alpha1_c,upstroke,downstroke_beginning,theta,theta_min,R,RD,S,P,T,fprime_n_i,fprime_m_i,fprime_c_i,fprime_n_db_i,fprime_m_db_i,fprime_c_db_i,T_db,T_s,RD_tv0,T_flag_downstroke,alpha1_0n,alpha1_0m,alpha1_0c,beta_Sig1n,beta_Sig1c,beta_Sig2n,beta_S2n_lpr,beta_S2c_lpr,beta_S1n_u,beta_S1m_u,beta_S1c_u,beta_S1n_d,beta_S1m_d,beta_S1c_d,beta_S2n_u,beta_S2m_u,beta_S2c_u,beta_S2n_d,beta_S2m_d,beta_S2c_d,xi,f0_n,f0_m,f0_c,fb_n,fb_m,fb_c,S1_n,S1_m,S1_c,S2_n,S2_m,S2_c);

%% Find time of stall onset
[tv0,tau_v,f_diff_tv0,qR_tv0,R_tv0,RD_tv0,Tv_tv0,theta_tv0,upstroke_tv0,f_diff_tv0_2,RD_tv0_2,upstroke_tv0_2] = BL_stall_time(tv0,f_diff_tv0,qR_tv0,R_tv0,RD_tv0,Tv_tv0,theta_tv0,upstroke_tv0,f_diff_tv0_2,RD_tv0_2,upstroke_tv0_2,t_i,theta_i,t,theta,upstroke,qR,R,RD,f_n,f2prime_n,g_v,Tv);

%% Time delay variables
[Tf_n,Tf_m,Tf_c,Ta_theta] = BL_time_constants(in_stall,upstroke,theta,qR,R,RD,P,qR_max,lambda_1,lambda_2,Ta,Tf);

%% Airloads coefficients
% Vortex overshoots  
[c_nv,c_mv,c_cv] = BL_vortex_overshoots(tau_v,Tv_tv0,f_diff_tv0,qR_tv0,R_tv0,RD_tv0,theta_tv0,upstroke_tv0,f_diff_tv0_2,RD_tv0_2,upstroke_tv0_2,mu_v2,nu_1,nu_2,nu_3,nu_4,nu_5,chi_u,chi_d,g_v,g_v2,Tv,Tv2,Vm,Vc,Vn1,Vn2,Vn3);
% c_n
[c_n,c_nC,c_nf,c_nI,c_nIa,c_nIf,c_nIg,w_E,w_Ep,w_Ef,w_Eg,alpha_E] = BLhargen_cn_coeff(x,gust_states,c_nv,f2prime_n,U,b,beta,Uc,w_tqc_pp,deltadot,deltaddot,w_tqc_f,wg,vg,AG,bG,G,Cg,eps_fn,alpha_0L,c_n_alpha,T1,T4);
% c_m
[c_m,c_mC,c_mf,c_mI,c_mQS,dCP] = BLhargen_cm_coeff(airfoil,x,c_nC,c_nf,c_mv,Ucm_QS,f2prime_m,theta,R,S,P,RD,upstroke,ah,dh,U,b,M,delta,deltadot,deltaddot,eps_fm,kappa_0,kappa_1,kappa_2,kappa_3,c_m0,T1,T4,T7,T8,T10,T11);
% c_c
c_c = BL_cc_coeff(c_cv,upstroke,alpha,alpha_lag,theta,f2prime_c,alpha_E,R,RD,S,T_s,alpha_0L,alpha1_0c,eta,c_d0,c_n_alpha,E0,E1);
% c_h
if any([d_0,d_1])
    c_h = BLTgen_ch_coeff(U,b,ah,w_Ep,w_Ef,alphadot,alphaddot,hddot,delta,deltadot,deltaddot,T1,T3,T4,T5,T9,T10,T11,T12,T13,alpha_0L,eps_fh);
else
    c_h = 0;
end

%% Outputs
y = [alpha, alpha_plunge, delta, w_QS, alpha_QS, q_QS, qR, R, alpha_cr, theta, theta_min, theta_max, S, P, T, alpha1_n, alpha1_m, alpha1_c, dalpha1_n, dalpha1_m, dalpha1_c, f_n, f_m, f_c, fprime_n, fprime_m, fprime_c, Tf_n, Tf_m, Tf_c, Ta_theta, alpha_E, c_n, c_nC, c_nI, c_nf, c_nv, c_m, c_mC, c_mI, c_mf, c_mv, dCP, c_c, c_h, RD_tv0, f_diff_tv0, Tv_tv0, f_diff_tv0_2, qR_max, U, beta, Uc, vg, wg, w_tqc_pp, wdot_tqc_pp, rdot, T_s, w_E, c_cv, c_mQS];

end