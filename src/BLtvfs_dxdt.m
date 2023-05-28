function [xdot,tv0,t_db,alpha_lag,so_lim,alpha1_n,alpha1_m,alpha1_c,alpha,q,fprime_n,fprime_m,fprime_c,dalpha1_n,beta_db_n,beta_db_m,beta_db_c,dalpha1_db,RD_db,R,g_v_tv0,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F] = ...
          BLtvfs_dxdt(t,x,xdot,t_i,tv0,t_db,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F,so_i,so_lim_i,upstroke_i,fprime_n_i,fprime_m_i,fprime_c_i,dalpha1_n_i,beta_db_n,beta_db_m,beta_db_c,dalpha1_db,RD_db,a_inf,b,a_0,a_1,k,U_0,U_1,k_U,psi_a,A1,A2,A3,A4,b1,b2,b3,b4,b5,S1,S2,TvL,Ta,Tf0,r0,alpha_ds0,alpha_ss,alpha1_0,f0,fb,delta_alpha_0,delta_alpha_1,delta_alpha_2,d_cm,d_cc,z_cm,z_cc,gamma_LS,g_v,df0_c,fSS1,fSS2,fSig1n,fSig1c,fSig2n,fS2n_ulpr,fS2n_dlpr,fS1n_u,fS1m_u,fS1c_u,fS1n_d,fS1m_d,fS1c_d,fS1n_ud,fS1m_ud,fS1c_ud,fS2n_u,fS2m_u,fS2c_u,fS2n_d,fS2m_d,fS2c_d,gamma_TvL,theta_max,theta_min,RD_m,R_max)

%% Kinematics
[U,M,Mdot,beta,alpha,alphadot,alphaddot,w_tqc,w_tqc_dot,cnaw_tqc_dot,q,qdot,r,rdot,qR,qRdot,R,Rdot,Ucm_qs,Ucmqs_dot] = BLtvfs_kinematics(t,U_0,U_1,k_U,b,k,a_0,a_1,r0,a_inf,psi_a);

%% States
alpha_lag = x(1);
f2prime_n = x(2);
RD = x(5);
RD_theta = x(6);

%% Stall onset criterion
[so_lim,theta,S,P,T,upstroke,downstroke_beginning,T_db,t_db,RD_db] = BL_stall_onset(t,t_db,RD_db,upstroke_i,q,theta_max,alpha_lag,RD,RD_theta,RD_m,alpha_ds0,alpha_ss,gamma_LS,Tf0);

%% Unsteady breakpoint of separation angles, stall and motion qualifiers
[alpha1_n,dalpha1_n,alpha1_m,dalpha1_m,alpha1_c,dalpha1_c,dalpha1_db] = BL_alpha_brk(T_db,upstroke,rdot,R_max,qR,qRdot,R,RD,RD_theta,S,P,T,alpha1_0,alpha_ds0,alpha_ss,delta_alpha_0,delta_alpha_1,delta_alpha_2,d_cm,d_cc,z_cm,z_cc,downstroke_beginning,dalpha1_n_i,dalpha1_db);
    
%% Separation points
[f,fprime_n,fprime_m,fprime_c,Sigma2,beta_db_n,beta_db_m,beta_db_c] = BL_sep_points(alpha,alpha_lag,alpha1_n,alpha1_m,alpha1_c,upstroke,downstroke_beginning,theta,theta_min,R,RD,Rdot,S,P,T,fprime_n_i,fprime_m_i,fprime_c_i,beta_db_n,beta_db_m,beta_db_c,T_db,RD_tv0,alpha1_0,S1,S2,f0,fb,df0_c,fSig1n,fSig1c,fSig2n,fS2n_ulpr,fS2n_dlpr,fS1n_u,fS1m_u,fS1c_u,fS1n_d,fS1m_d,fS1c_d,fS1n_ud,fS1m_ud,fS1c_ud,fS2n_u,fS2m_u,fS2c_u,fS2n_d,fS2m_d,fS2c_d);

%% Find time of stall onset
[tv0,~,g_v_tv0,f_diff_tv0,RD_tv0,TvL_tv0,theta_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F] = BL_stall_time(tv0,f_diff_tv0,RD_tv0,TvL_tv0,theta_tv0,alpha_lag,so_lim,so_i,so_lim_i,t,t_i,f2prime_n,f,TvL,R,RD,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F,g_v,theta,gamma_TvL,q);

%% Time delay constants variation
[Tf_n,Tf_m,Tf_c,Ta_theta] = BL_time_constants(upstroke,theta,R,RD,P,S,Ta,Tf0,fSS1,fSS2);

%% State-space
xdot = BLtvfs_state_space(x,xdot,U,M,b,beta,Mdot,alpha,alphadot,alphaddot,cnaw_tqc_dot,Ucmqs_dot,R,Tf_n,Tf_m,Tf_c,Ta,fprime_n,fprime_m,fprime_c,Ta_theta,A1,A2,A3,A4,b1,b2,b3,b4,b5);

end