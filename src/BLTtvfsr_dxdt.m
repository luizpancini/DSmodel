function [xdot,tv0,alpha_lag,so_lim,alpha1_n,alpha1_m,alpha1_c,alpha,q,g_v_tv0,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F] = ...
         BLTtvfsr_dxdt(t,x,xdot,t_i,tv0,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F,so_i,so_lim_i,b1,b2,a_inf,U_0,U_1,k_U,b,a_0,a_1,k,psi_a,S1,S2,TvL,Ta,Tf0,r0,alpha_ds0,alpha_ss,alpha1_0,f0,fb,delta_alpha_0,delta_alpha_1,delta_alpha_2,d_cm,d_cc,z_cm,z_cc,gamma_LS,g_v,df0_c,fSS1,fSS2,fSig1n,fSig1c,fSig2n,fS2n_ulpr,fS2n_dlpr,fS1n_u,fS1m_u,fS1c_u,fS1n_d,fS1m_d,fS1c_d,fS1n_ud,fS1m_ud,fS1c_ud,fS2n_u,fS2m_u,fS2c_u,fS2n_d,fS2m_d,fS2c_d,gamma_TvL,theta_max,theta_min,RD_m)

%% Kinematics
[~,~,~,~,alpha,~,~,~,~,q,~,~,~,qR,qRdot,R,Rdot] = tvfs_kinematics(t,U_0,U_1,k_U,b,k,a_0,a_1,r0,a_inf,psi_a);

%% States
alpha_lag = x(1);
f2prime_n = x(2);
RD = x(5);
RD_theta = x(6);

%% Stall onset criterion
[so_lim,theta] = BL_stall_onset(alpha_lag,alpha_ss,alpha_ds0,RD_theta);         

%% Unsteady breakpoint of separation angles and stall qualifiers
[alpha1_n,~,alpha1_m,~,alpha1_c,~,P,S,T] = BL_alpha_brk(alpha1_0,alpha_ss,alpha_ds0,delta_alpha_0,delta_alpha_1,delta_alpha_2,d_cm,d_cc,z_cm,z_cc,R,RD,qR,q,theta,theta_max,RD_m,RD_theta,gamma_LS,qRdot);

%% Separation points
[f,fprime_n,fprime_m,fprime_c] = BL_sep_points(alpha_lag,f0,fb,alpha,alpha1_n,alpha1_0,S1,S2,R,q,theta,RD,alpha1_m,alpha1_c,theta_max,theta_min,RD_tv0,S,P,T,df0_c,fSig1n,fSig1c,fSig2n,fS2n_ulpr,fS2n_dlpr,fS1n_u,fS1m_u,fS1c_u,fS1n_d,fS1m_d,fS1c_d,fS1n_ud,fS1m_ud,fS1c_ud,fS2n_u,fS2m_u,fS2c_u,fS2n_d,fS2m_d,fS2c_d);

%% Find time of stall onset
[tv0,~,g_v_tv0,f_diff_tv0,RD_tv0,TvL_tv0,theta_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F] = BL_stall_time(tv0,f_diff_tv0,RD_tv0,TvL_tv0,theta_tv0,alpha_lag,so_lim,so_i,so_lim_i,t,t_i,f2prime_n,f,TvL,R,RD,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F,g_v,theta,gamma_TvL,q);

%% Time delay constants variation
[Tf_n,Tf_m,Tf_c,Ta_theta] = BL_time_constants(Ta,Tf0,theta,R,q,P,RD_theta,S,RD,fSS1,fSS2);

%% State-space
xdot = BLTtvfsr_state_space(x,xdot,alpha,R,Tf_n,Tf_m,Tf_c,Ta,fprime_n,fprime_m,fprime_c,Ta_theta);

end