function [y,F_theta0,F2_theta0] = BLTgust_outputs(t,dt,x,F_theta0_p,F2_theta0_p,tv0,U,b,beta,k,a_0,a_1,A1,A2,b1,b2,Cg,Cg2,lambda_g,AG,K0,K1,K1_f,K2,K3,kappa,alpha_0L,E0,E1,c_m0,c_n_alpha,r0,alpha_ds0,alpha_ss,alpha1_0,S1,S2,Ta,Tf0,f0,fb,Vn1,Vn2,Vm,c_d0,delta_alpha_0,delta_alpha_1,delta_alpha_2,d_cm,d_cc,z_cm,z_cc,nu_1,nu_2,gamma_LS,df0_c,fSig1n,fSig1c,fSig2n,fS2n_ulpr,fS2n_dlpr,fS1n_u,fS1m_u,fS1c_u,fS1n_d,fS1m_d,fS1c_d,fS1n_ud,fS1m_ud,fS1c_ud,fS2n_u,fS2m_u,fS2c_u,fS2n_d,fS2m_d,fS2c_d,fSS1,fSS2,g_v_tv0,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,theta_max,theta_min,RD_m,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,TvL,wg_fun,wgdot_fun)

%% Kinematics
[alpha_g,alpha_g0,alpha,alphadot,alphaddot,alpha_tqc,q,~,~,~,qR,qRdot,R,Rdot,F_theta0,F2_theta0] = convecting_gust_kinematics(t,U,b,k,a_0,a_1,r0,wg_fun,wgdot_fun,lambda_g);

%% States
alpha_hat_lag = x(3)+x(end); % Composition of motion-induced and gust-induced AoA lagged
f2prime_n = x(4);
f2prime_m = x(5);
f2prime_c = x(6);
RD = x(7);
RD_theta = x(8);
gust_states = x(9:end-1);

%% Stall onset criterion
[so_lim,theta] = BL_stall_onset(alpha_hat_lag,alpha_ss,alpha_ds0,RD_theta);         

%% Unsteady breakpoint of separation angles and stall qualifiers
[alpha1_n,dalpha1_n,alpha1_m,dalpha1_m,alpha1_c,dalpha1_c,P,S,T] = BL_alpha_brk(alpha1_0,alpha_ss,alpha_ds0,delta_alpha_0,delta_alpha_1,delta_alpha_2,d_cm,d_cc,z_cm,z_cc,R,RD,qR,q,theta,theta_max,RD_m,RD_theta,gamma_LS,qRdot);

%% Separation points
[f,fprime_n,fprime_m,fprime_c,Sigma2] = BL_sep_points(alpha_hat_lag,f0,fb,alpha,alpha1_n,alpha1_0,S1,S2,R,q,theta,RD,alpha1_m,alpha1_c,theta_max,theta_min,RD_tv0,S,P,T,df0_c,fSig1n,fSig1c,fSig2n,fS2n_ulpr,fS2n_dlpr,fS1n_u,fS1m_u,fS1c_u,fS1n_d,fS1m_d,fS1c_d,fS1n_ud,fS1m_ud,fS1c_ud,fS2n_u,fS2m_u,fS2c_u,fS2n_d,fS2m_d,fS2c_d);

%% Find time of stall onset
tau_v = max([0 t-tv0]);

%% Time delay constants variation
[Tf_n,Tf_m,Tf_c] = BL_time_constants(Ta,Tf0,theta,R,q,P,RD_theta,S,RD,fSS1,fSS2);

%% Airloads coefficients
% Vortex overshoots  
[cn_v,cm_v,cc_v] = BL_vortex_overshoots(tau_v,TvL_tv0,f_diff_tv0,RD_tv0,theta_tv0,Vn1,Vn2,Vm,nu_1,nu_2,g_v_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,q,P,R,Sigma2,theta,TvL);
% c_n
[c_n,c_nC,c_nf,c_nI,alpha_Et,K_f] = BLTgust_cn_coeff(x,gust_states,U,b,beta,alphadot,alphaddot,alpha_tqc,A1,A2,b1,b2,Cg,Cg2,AG,c_n_alpha,f2prime_n,alpha_0L,cn_v,alpha_g,alpha_g0,F_theta0,F_theta0_p,dt);
% c_m
[c_m,c_mf,c_mI,dCP] = BLTgust_cm_coeff(c_nf,alphadot,alphaddot,q,K0,K1,K1_f,K2,K3,kappa,c_m0,U,b,f2prime_m,theta,R,S,P,RD,cm_v,alpha_g0,F2_theta0,F2_theta0_p,dt);
% c_c
c_c = BL_cc_coeff(alpha,f2prime_c,theta,c_n_alpha,alpha_Et,E0,E1,c_d0,alpha_0L,RD,cc_v);
 
%% Outputs
y = [alpha; c_n; c_m; c_c; c_nf; c_nI; c_mf; c_mI; cn_v; cm_v; cc_v; f; tau_v; dCP; so_lim; qR; alpha1_n; K_f; c_nC; Tf_n; dalpha1_n; fprime_n; fprime_m; fprime_c; q; dalpha1_m; dalpha1_c; Rdot; Tf_m; Tf_c; theta_min; theta_max; P; S; alpha_Et; alpha_g; U];

end