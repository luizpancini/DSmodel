function [y,cnaw_tqc,alpha,alphadot,M,Ucm_qs,s,q] = BLtvfsc_outputs(t,x,tv0,tvec,cnaw_tqc_vec,alpha_vec,alphadot_vec,M_vec,Ucm_qs_vec,s_vec,qdot_vec,a_inf,b,k,a_0,a_1,U_0,U_1,k_U,psi_a,K0,K1,K1_f,K2,K3,kappa,b1,b2,b3,b4,b5,A1,A2,A3,A4,alpha_0L,E0,E1,c_m0,r0,alpha_ds0,alpha_ss,alpha1_0,S1,S2,Ta,Tf0,f0,fb,Vn1,Vn2,Vm,c_d0,delta_alpha_0,delta_alpha_1,delta_alpha_2,d_cm,d_cc,z_cm,z_cc,nu_1,nu_2,gamma_LS,df0_c,fSig1n,fSig1c,fSig2n,fS2n_ulpr,fS2n_dlpr,fS1n_u,fS1m_u,fS1c_u,fS1n_d,fS1m_d,fS1c_d,fS1n_ud,fS1m_ud,fS1c_ud,fS2n_u,fS2m_u,fS2c_u,fS2n_d,fS2m_d,fS2c_d,fSS1,fSS2,g_v_tv0,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,theta_max,theta_min,RD_m,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,TvL)

%% Kinematics
[U,M,~,beta,alpha,alphadot,~,w_tqc,~,~,q,~,~,~,qR,qRdot,R,Rdot,Ucm_qs,~] = BLtvfs_kinematics(t,U_0,U_1,k_U,b,k,a_0,a_1,r0,a_inf,psi_a);

%% States
alpha_lag = x(1);
f2prime_n = x(2);
f2prime_m = x(3);
f2prime_c = x(4);
RD = x(5);
RD_theta = x(6);

%% Stall onset criterion
[so_lim,theta] = BL_stall_onset(alpha_lag,alpha_ss,alpha_ds0,RD_theta);         

%% Unsteady breakpoint of separation angles and stall qualifiers
[alpha1_n,dalpha1_n,alpha1_m,dalpha1_m,alpha1_c,dalpha1_c,P,S,T] = BL_alpha_brk(alpha1_0,alpha_ss,alpha_ds0,delta_alpha_0,delta_alpha_1,delta_alpha_2,d_cm,d_cc,z_cm,z_cc,R,RD,qR,q,theta,theta_max,RD_m,RD_theta,gamma_LS,qRdot);

%% Separation points
[f,fprime_n,fprime_m,fprime_c,Sigma2] = BL_sep_points(alpha_lag,f0,fb,alpha,alpha1_n,alpha1_0,S1,S2,R,q,theta,RD,alpha1_m,alpha1_c,theta_max,theta_min,RD_tv0,S,P,T,df0_c,fSig1n,fSig1c,fSig2n,fS2n_ulpr,fS2n_dlpr,fS1n_u,fS1m_u,fS1c_u,fS1n_d,fS1m_d,fS1c_d,fS1n_ud,fS1m_ud,fS1c_ud,fS2n_u,fS2m_u,fS2c_u,fS2n_d,fS2m_d,fS2c_d);

%% Find time of stall onset
tau_v = max([0 t-tv0]);

%% Time delay constants variation
[Tf_n,Tf_m,Tf_c] = BL_time_constants(Ta,Tf0,theta,R,q,P,RD_theta,S,RD,fSS1,fSS2);

%% Airloads coefficients
% Airspeed-dependent parameters
[c_n_alpha,K0] = airspeed_vars(M,beta);
% Convolution integrals
[int_conv_Ia,int_conv_Iq,int_conv_IM,int_conv_Iam,int_conv_Iqm,int_conv_IMm,cnaw_tqc,w_tqc_eff,Ucm_qs_eff,s] = BLtvfsc_convolution_integrals(t,a_inf,b,U,M,beta,alpha,alphadot,q,w_tqc,Ucm_qs,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,b5,tvec,cnaw_tqc_vec,alpha_vec,alphadot_vec,M_vec,Ucm_qs_vec,s_vec,qdot_vec);
% Vortex overshoots  
[cn_v,cm_v,cc_v] = BL_vortex_overshoots(tau_v,TvL_tv0,f_diff_tv0,RD_tv0,theta_tv0,Vn1,Vn2,Vm,nu_1,nu_2,g_v_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,q,P,R,Sigma2,theta,TvL);
% c_n
[c_n,c_nC,c_nf,c_nI,alpha_E,K_f] = BLtvfsc_cn_coeff(U,w_tqc_eff,int_conv_Ia,int_conv_Iq,int_conv_IM,c_n_alpha,f2prime_n,alpha_0L,cn_v);
% c_m
[c_m,c_mf,c_mI,dCP] = BLtvfsc_cm_coeff(c_nf,U,q,Ucm_qs_eff,int_conv_Iam,int_conv_Iqm,int_conv_IMm,K0,K1,K1_f,K2,K3,kappa,c_m0,f2prime_m,theta,R,S,P,RD,cm_v);
% c_c
c_c = BL_cc_coeff(alpha,f2prime_c,theta,c_n_alpha,alpha_E,E0,E1,c_d0,alpha_0L,RD,cc_v);
 
%% Outputs
y = [alpha; c_n; c_m; c_c; c_nf; c_nI; c_mf; c_mI; cn_v; cm_v; cc_v; f; tau_v; dCP; so_lim; qR; alpha1_n; K_f; c_nC; Tf_n; dalpha1_n; fprime_n; fprime_m; fprime_c; q; dalpha1_m; dalpha1_c; Rdot; Tf_m; Tf_c; theta_min; theta_max; P; S; alpha_E; U; beta];

end