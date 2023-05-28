function y = BLtvfs_outputs(t,x,tv0,t_db,upstroke_i,fprime_n_i,fprime_m_i,fprime_c_i,dalpha1_n_i,beta_db_n,beta_db_m,beta_db_c,dalpha1_db,RD_db,a_inf,b,k,a_0,a_1,U_0,U_1,k_U,psi_a,K0,K1,K1_f,K2,K3,kappa,alpha_0L,E0,E1,c_m0,r0,alpha_ds0,alpha_ss,alpha1_0,S1,S2,Ta,Tf0,f0,fb,Vn1,Vn2,Vm,c_d0,delta_alpha_0,delta_alpha_1,delta_alpha_2,d_cm,d_cc,z_cm,z_cc,nu_1,nu_2,gamma_LS,df0_c,fSig1n,fSig1c,fSig2n,fS2n_ulpr,fS2n_dlpr,fS1n_u,fS1m_u,fS1c_u,fS1n_d,fS1m_d,fS1c_d,fS1n_ud,fS1m_ud,fS1c_ud,fS2n_u,fS2m_u,fS2c_u,fS2n_d,fS2m_d,fS2c_d,fSS1,fSS2,g_v_tv0,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,theta_max,theta_min,RD_m,R_max,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,TvL)

%% Kinematics
[U,M,Mdot,beta,alpha,alphadot,alphaddot,w_tqc,w_tqc_dot,cnaw_tqc_dot,q,qdot,r,rdot,qR,qRdot,R,Rdot,Ucm_qs,Ucmqs_dot] = BLtvfs_kinematics(t,U_0,U_1,k_U,b,k,a_0,a_1,r0,a_inf,psi_a);

%% States
alpha_lag = x(1);
f2prime_n = x(2);
f2prime_m = x(3);
f2prime_c = x(4);
RD = x(5);
RD_theta = x(6);

%% Stall onset criterion
[so_lim,theta,S,P,T,upstroke,downstroke_beginning,T_db] = BL_stall_onset(t,t_db,RD_db,upstroke_i,q,theta_max,alpha_lag,RD,RD_theta,RD_m,alpha_ds0,alpha_ss,gamma_LS,Tf0);

%% Unsteady breakpoint of separation angles, stall and motion qualifiers
[alpha1_n,dalpha1_n,alpha1_m,dalpha1_m,alpha1_c,dalpha1_c] = BL_alpha_brk(T_db,upstroke,rdot,R_max,qR,qRdot,R,RD,RD_theta,S,P,T,alpha1_0,alpha_ds0,alpha_ss,delta_alpha_0,delta_alpha_1,delta_alpha_2,d_cm,d_cc,z_cm,z_cc,downstroke_beginning,dalpha1_n_i,dalpha1_db);

%% Separation points
[f,fprime_n,fprime_m,fprime_c,Sigma2] = BL_sep_points(alpha,alpha_lag,alpha1_n,alpha1_m,alpha1_c,upstroke,downstroke_beginning,theta,theta_min,R,RD,Rdot,S,P,T,fprime_n_i,fprime_m_i,fprime_c_i,beta_db_n,beta_db_m,beta_db_c,T_db,RD_tv0,alpha1_0,S1,S2,f0,fb,df0_c,fSig1n,fSig1c,fSig2n,fS2n_ulpr,fS2n_dlpr,fS1n_u,fS1m_u,fS1c_u,fS1n_d,fS1m_d,fS1c_d,fS1n_ud,fS1m_ud,fS1c_ud,fS2n_u,fS2m_u,fS2c_u,fS2n_d,fS2m_d,fS2c_d);

%% Find time of stall onset
tau_v = max([0 t-tv0]);

%% Time delay constants variation
[Tf_n,Tf_m,Tf_c,Ta_theta] = BL_time_constants(upstroke,theta,R,RD,P,S,Ta,Tf0,fSS1,fSS2);

%% Airloads coefficients
% Airspeed-dependent parameters
[c_n_alpha,K0] = airspeed_vars(M,beta);
% Vortex overshoots  
[cn_v,cm_v,cc_v] = BL_vortex_overshoots(tau_v,TvL_tv0,f_diff_tv0,RD_tv0,theta_tv0,Vn1,Vn2,Vm,nu_1,nu_2,g_v_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,q,P,R,Sigma2,theta,TvL);
% c_n
[c_n,c_nC,c_nf,c_nI,alpha_E,K_f] = BLtvfs_cn_coeff(x,U,w_tqc,c_n_alpha,f2prime_n,alpha_0L,cn_v);
% c_m
[c_m,c_mf,c_mI,dCP] = BLtvfs_cm_coeff(x,c_nf,U,q,Ucm_qs,K0,K1,K1_f,K2,K3,kappa,c_m0,f2prime_m,theta,R,S,P,RD,cm_v,alphadot,c_n_alpha,b);
% c_c
c_c = BL_cc_coeff(alpha,f2prime_c,theta,c_n_alpha,alpha_E,E0,E1,c_d0,alpha_0L,RD,cc_v);
 
%% Outputs
y = [alpha; c_n; c_m; c_c; c_nf; c_nI; c_mf; c_mI; cn_v; cm_v; cc_v; f; tau_v; dCP; so_lim; qR; alpha1_n; K_f; c_nC; Tf_n; dalpha1_n; fprime_n; fprime_m; fprime_c; q; dalpha1_m; dalpha1_c; Rdot; Tf_m; Tf_c; theta_min; theta_max; P; S; alpha_E; U; beta];

end