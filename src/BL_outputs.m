function [y,t_ub,t_db,qR,theta,alpha_lag,dalpha1_n,fprime_n,fprime_m,fprime_c,fprime_n_db,fprime_m_db,fprime_c_db,RD_ub,RD_db,theta_min,theta_max,RD_m,qR_max,S,P,upstroke,T_flag_downstroke,tv0,f_diff_tv0,qR_tv0,R_tv0,RD_tv0,Tv_tv0,theta_tv0,upstroke_tv0,f_diff_tv0_2,RD_tv0_2,upstroke_tv0_2,r] = ...
         BL_outputs(t,x,tv0,t_ub,t_db,t_i,dt,r_i,qR_i,theta_i,alpha_lag_i,upstroke_i,S_i,P_i,T_flag_downstroke_i,dalpha1_n_i,dalpha1_db_i,fprime_n_i,fprime_m_i,fprime_c_i,fprime_n_db_i,fprime_m_db_i,fprime_c_db_i,RD_ub_i,RD_db_i,theta_max_i,theta_min_i,RD_m_i,qR_max_i,qR_tv0,R_tv0,RD_tv0,f_diff_tv0,theta_tv0,upstroke_tv0,Tv_tv0,f_diff_tv0_2,RD_tv0_2,upstroke_tv0_2,b,ah,U,M,beta,k,a_0,a_1,k_h,a_1h,K_a,K_aM,K_q,K_qM,T_I,b1,b2,b3,b4,b5,A1,A2,A3,A4,alpha_0L,alpha1_0n,alpha1_0m,alpha1_0c,alpha_ds0,alpha_ss,beta_Sig1n,beta_Sig1c,beta_Sig2n,beta_S2n_lpr,beta_S2c_lpr,beta_S1n_u,beta_S1m_u,beta_S1c_u,beta_S1n_d,beta_S1m_d,beta_S1c_d,beta_S2n_u,beta_S2m_u,beta_S2c_u,beta_S2n_d,beta_S2m_d,beta_S2c_d,gamma_LS,delta_alpha_0,delta_alpha_1,eta,kappa_0,kappa_1,kappa_2,kappa_3,lambda_1,lambda_2,mu_v2,nu_1,nu_2,nu_3,nu_4,nu_5,chi_u,chi_d,xi,zeta_a,c_d0,c_m0,c_n_alpha,d_cc,d_cm,E0,E1,f0_n,f0_m,f0_c,fb_n,fb_m,fb_c,g_v,g_v2,K0,K1,K2,r0,S1_n,S1_m,S1_c,S2_n,S2_m,S2_c,Ta,Tf,Tv,Tv2,Vn1,Vn2,Vn3,Vm,Vc,z_ccd,z_ccu,z_cm)

%% Kinematics
[alpha,alpha_plunge,w_QS,Uc,alpha_QS,hdot,w_tqc,q,q_QS,qR,R,r,rdot] = BL_kinematics(t,U,b,ah,k,a_0,a_1,a_1h,k_h,r0,r_i,dt);

%% States
alpha_lag = atan(x(9)/Uc);
f2prime_n = x(10);
f2prime_m = x(11);
f2prime_c = x(12);
RD = x(13);
RD_theta = x(14);

%% Stall onset criterion, stall and motion qualifiers
[alpha_cr,theta,theta_max,theta_min,del_RD_acc,RD_m,qR_max,R_max,S,P,T,upstroke,downstroke_beginning,in_stall,T_flag_downstroke,t_ub,t_db,T_ub,T_db,RD_ub,RD_db,T_s] = BL_stall_onset(t,t_ub,t_db,tv0,RD_ub_i,RD_db_i,RD_tv0,theta_i,alpha_lag_i,qR_i,upstroke_i,S_i,P_i,T_flag_downstroke_i,theta_max_i,theta_min_i,RD_m_i,qR_max_i,alpha_lag,q_QS,qR,R,RD,RD_theta,rdot,alpha_ds0,alpha_ss,gamma_LS,Tf,zeta_a);

%% Unsteady breakpoint of separation angles
[alpha1_n,dalpha1_n,alpha1_m,dalpha1_m,alpha1_c,dalpha1_c] = BL_alpha_brk(qR,R,RD,RD_theta,del_RD_acc,RD_m,theta,theta_min,theta_max,R_max,S,P,T,T_s,T_ub,T_db,upstroke,downstroke_beginning,dalpha1_n_i,dalpha1_db_i,alpha1_0n,alpha1_0m,alpha1_0c,alpha_ds0,alpha_ss,delta_alpha_0,delta_alpha_1,d_cm,d_cc,z_cm,z_ccd,z_ccu);

%% Separation points
[f_n,f_m,f_c,fprime_n,fprime_m,fprime_c,Sigma2,fprime_n_db,fprime_m_db,fprime_c_db] = BL_sep_points(alpha_QS,alpha_lag,alpha1_n,alpha1_m,alpha1_c,upstroke,downstroke_beginning,theta,theta_min,R,RD,S,P,T,fprime_n_i,fprime_m_i,fprime_c_i,fprime_n_db_i,fprime_m_db_i,fprime_c_db_i,T_db,T_s,RD_tv0,T_flag_downstroke,alpha1_0n,alpha1_0m,alpha1_0c,beta_Sig1n,beta_Sig1c,beta_Sig2n,beta_S2n_lpr,beta_S2c_lpr,beta_S1n_u,beta_S1m_u,beta_S1c_u,beta_S1n_d,beta_S1m_d,beta_S1c_d,beta_S2n_u,beta_S2m_u,beta_S2c_u,beta_S2n_d,beta_S2m_d,beta_S2c_d,xi,f0_n,f0_m,f0_c,fb_n,fb_m,fb_c,S1_n,S1_m,S1_c,S2_n,S2_m,S2_c);

%% Find time of stall onset
[tv0,tau_v,f_diff_tv0,qR_tv0,R_tv0,RD_tv0,Tv_tv0,theta_tv0,upstroke_tv0,f_diff_tv0_2,RD_tv0_2,upstroke_tv0_2] = BL_stall_time(tv0,f_diff_tv0,qR_tv0,R_tv0,RD_tv0,Tv_tv0,theta_tv0,upstroke_tv0,f_diff_tv0_2,RD_tv0_2,upstroke_tv0_2,t_i,theta_i,t,theta,upstroke,qR,R,RD,f_n,f2prime_n,g_v,Tv);

%% Time delay variables
[Tf_n,Tf_m,Tf_c,Ta_theta] = BL_time_constants(in_stall,upstroke,theta,qR,R,RD,P,qR_max,lambda_1,lambda_2,Ta,Tf);

%% Airloads coefficients
% Vortex overshoots  
[c_nv,c_mv,c_cv] = BL_vortex_overshoots(tau_v,Tv_tv0,f_diff_tv0,qR_tv0,R_tv0,RD_tv0,theta_tv0,upstroke_tv0,f_diff_tv0_2,RD_tv0_2,upstroke_tv0_2,mu_v2,nu_1,nu_2,nu_3,nu_4,nu_5,chi_u,chi_d,g_v,g_v2,Tv,Tv2,Vm,Vc,Vn1,Vn2,Vn3);
% c_n
[c_n,c_nC,c_nf,c_nI,alpha_C] = BL_cn_coeff(x,M,U,b,ah,beta,Uc,alpha_QS,q,alpha,w_tqc,hdot,f2prime_n,alpha_0L,c_n_alpha,A1,A2,b1,b2,K_a,K_q,T_I,c_nv);
% c_m
[c_m,c_mC,c_mf,c_mI,dCP] = BL_cm_coeff(c_mv,c_nC,c_nf,x,f2prime_m,upstroke,M,U,b,ah,beta,alpha_QS,q,theta,R,RD,S,P,kappa_0,kappa_1,kappa_2,kappa_3,c_m0,K0,K1,K2,A3,b3,A4,b4,b5,K_aM,K_qM,T_I);
% c_c
c_c = BL_cc_coeff(c_cv,upstroke,alpha,alpha_lag,theta,f2prime_c,alpha_C,R,RD,S,T_s,alpha_0L,alpha1_0c,eta,c_d0,c_n_alpha,E0,E1);
 
%% Outputs
y = [alpha, alpha_plunge, alpha_QS, q_QS, qR, R, alpha_cr, theta, theta_min, theta_max, S, P, T, alpha1_n, alpha1_m, alpha1_c, dalpha1_n, dalpha1_m, dalpha1_c, f_n, f_m, f_c, fprime_n, fprime_m, fprime_c, Tf_n, Tf_m, Tf_c, Ta_theta, alpha_C, c_n, c_nC, c_nI, c_nf, c_nv, c_m, c_mC, c_mI, c_mf, c_mv, dCP, c_c, RD_tv0, f_diff_tv0, Tv_tv0, f_diff_tv0_2, qR_max, Uc];

end