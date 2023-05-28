function y = BLT_outputs(t,x,tv0,t_ub,t_db,theta_i,alpha_lag_i,R_i,upstroke_i,S_i,T_flag_downstroke,dalpha1_n_i,dalpha1_db,fprime_n_i,fprime_m_i,fprime_c_i,fprime_n_db,fprime_m_db,fprime_c_db,RD_ub,RD_db,theta_max,theta_min,RD_m,R_max,U,b,ah,beta,k,a_0,a_1,a_1h,k_h,c1,c2,c3,c4,A1,A2,b1,b2,K0,K1,K1_f,K2,K3,kappa,alpha_0L,E0,E1,c_m0,c_n_alpha,r0,alpha_ds0,alpha_ss,alpha1_0,S1,S2,Ta,Tf0,f0,fb,Vn1,Vn2,Vm,c_d0,delta_alpha_0,delta_alpha_1,delta_alpha_2,d_cm,d_cc,z_cm,z_cc,nu_1,nu_2,nu_3,gamma_LS,df0_c,fSig1n,fSig1c,fSig2n,fS2n_ulpr,fS2n_dlpr,fS1n_u,fS1m_u,fS1c_u,fS1n_d,fS1m_d,fS1c_d,fS1n_ud,fS1m_ud,fS1c_ud,fS2n_u,fS2m_u,fS2c_u,fS2n_d,fS2m_d,fS2c_d,fSS1,fSS2,g_v_tv0,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,TvL)

%% Kinematics
[alpha,alphadot,alphaddot,hddot,alpha_plunge,alpha_bar,alpha_tqc,q,q_bar,qdot_bar,r,rdot,qR,qRdot,R,Rdot] = BL_kinematics(t,U,b,ah,k,a_0,a_1,a_1h,k_h,r0);

%% States
alpha_lag = x(3);
f2prime_n = x(4);
f2prime_m = x(5);
f2prime_c = x(6);
RD = x(7);
RD_theta = x(8);

%% Stall onset criterion, stall and motion qualifiers
[alpha_cr,theta,theta_max,theta_min,RD_m,R_max,S,P,T,upstroke,downstroke_beginning,in_stall,T_flag_downstroke,t_ub,t_db,T_ub,T_db,RD_ub,RD_db,T_s] = BL_stall_onset(t,t_ub,t_db,tv0,RD_ub,RD_db,RD_tv0,theta_i,alpha_lag_i,R_i,upstroke_i,S_i,T_flag_downstroke,theta_max,theta_min,RD_m,R_max,alpha_lag,q_bar,R,RD,RD_theta,alpha_ds0,alpha_ss,gamma_LS,Tf0);

%% Unsteady breakpoint of separation angles
[alpha1_n,dalpha1_n,alpha1_m,dalpha1_m,alpha1_c,dalpha1_c] = BL_alpha_brk(qR,qRdot,R,RD,RD_theta,theta,theta_max,R_max,S,P,T,T_s,T_ub,T_db,upstroke,downstroke_beginning,dalpha1_n_i,dalpha1_db,alpha1_0,alpha_ds0,alpha_ss,delta_alpha_0,delta_alpha_1,delta_alpha_2,d_cm,d_cc,z_cm,z_cc);

%% Separation points
[f,fprime_n,fprime_m,fprime_c,Sigma2] = BL_sep_points(alpha_bar,alpha_lag,alpha1_n,alpha1_m,alpha1_c,upstroke,downstroke_beginning,theta,theta_min,R,RD,Rdot,S,P,T,fprime_n_i,fprime_m_i,fprime_c_i,fprime_n_db,fprime_m_db,fprime_c_db,T_db,T_s,RD_tv0,T_flag_downstroke,alpha1_0,S1,S2,f0,fb,df0_c,fSig1n,fSig1c,fSig2n,fS2n_ulpr,fS2n_dlpr,fS1n_u,fS1m_u,fS1c_u,fS1n_d,fS1m_d,fS1c_d,fS1n_ud,fS1m_ud,fS1c_ud,fS2n_u,fS2m_u,fS2c_u,fS2n_d,fS2m_d,fS2c_d);

%% Find time of stall onset
tau_v = max([0 t-tv0]);

%% Time delay variables
[Tf_n,Tf_m,Tf_c,Ta_theta] = BL_time_constants(in_stall,upstroke,theta,R,RD,P,Ta,Tf0,fSS1,fSS2);

%% Airloads coefficients
% Vortex overshoots  
[c_nv,c_mv,c_cv] = BL_vortex_overshoots(tau_v,TvL_tv0,f_diff_tv0,RD_tv0,theta_tv0,Vn1,Vn2,Vm,nu_1,nu_2,nu_3,g_v_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,q_bar,P,R,Sigma2,theta,TvL);
% c_n
[c_n,c_nC,c_nf,c_nI,alpha_E,K_f] = BLT_cn_coeff(x,U,b,ah,beta,alphadot,alphaddot,hddot,alpha_tqc,c1,c2,c3,c4,A1,A2,b1,b2,c_n_alpha,f2prime_n,alpha_0L,c_nv);
% c_m
[c_m,c_mf,c_mI,dCP] = BLT_cm_coeff(c_nf,alphadot,alphaddot,hddot,q_bar,K0,K1,K1_f,K2,K3,kappa,c_m0,U,b,ah,f2prime_m,theta,R,S,P,RD,c_mv);
% c_c
c_c = BL_cc_coeff(in_stall,alpha,f2prime_c,theta,c_n_alpha,alpha_E,E0,E1,c_d0,alpha_0L,RD,P,c_cv);
 
%% Outputs
y = [alpha; c_n; c_m; c_c; c_nf; c_nI; c_mf; c_mI; c_nv; c_mv; c_cv; f; tau_v; dCP; alpha_cr; qR; alpha1_n; K_f; c_nC; Tf_n; dalpha1_n; fprime_n; fprime_m; fprime_c; q_bar; dalpha1_m; dalpha1_c; Rdot; Tf_m; Tf_c; theta_min; theta_max; P; S; alpha_E; alpha_plunge];

end