function [xdot,tv0,t_ub,t_db,alpha_bar,q_bar,R,alpha_lag,alpha_cr,alpha1_n,alpha1_m,alpha1_c,dalpha1_n,dalpha1_db,fprime_n,fprime_m,fprime_c,fprime_n_db,fprime_m_db,fprime_c_db,RD_ub,RD_db,g_v_tv0,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F] = ...
         BLT_dxdt(t,x,xdot,t_i,tv0,t_ub,t_db,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F,alpha_lag_i,alpha_cr_i,theta_i,R_i,upstroke_i,S_i,T_flag_downstroke_i,theta_max_i,theta_min_i,RD_m_i,R_max_i,dalpha1_n_i,dalpha1_db_i,fprime_n_i,fprime_m_i,fprime_c_i,fprime_n_db_i,fprime_m_db_i,fprime_c_db_i,RD_ub_i,RD_db_i,A,U,b,ah,a_0,a_1,k,a_1h,k_h,S1,S2,TvL,Ta,Tf0,r0,alpha_ds0,alpha_ss,alpha1_0,f0,fb,delta_alpha_0,delta_alpha_1,delta_alpha_2,d_cm,d_cc,z_cm,z_cc,gamma_LS,gamma_TvLs,gamma_TvLc,g_v,df0_c,fSS1,fSS2,fSig1n,fSig1c,fSig2n,fS2n_ulpr,fS2n_dlpr,fS1n_u,fS1m_u,fS1c_u,fS1n_d,fS1m_d,fS1c_d,fS1n_ud,fS1m_ud,fS1c_ud,fS2n_u,fS2m_u,fS2c_u,fS2n_d,fS2m_d,fS2c_d)

%% Kinematics
[alpha,alphadot,alphaddot,hddot,alpha_plunge,alpha_bar,alpha_tqc,q,q_bar,qdot_bar,r,rdot,qR,qRdot,R,Rdot] = BL_kinematics(t,U,b,ah,k,a_0,a_1,a_1h,k_h,r0);

%% States
alpha_lag = x(3);
f2prime_n = x(4);
RD = x(7);
RD_theta = x(8);

%% Stall onset criterion, stall and motion qualifiers
[alpha_cr,theta,theta_max,theta_min,RD_m,R_max,S,P,T,upstroke,downstroke_beginning,in_stall,T_flag_downstroke,t_ub,t_db,T_ub,T_db,RD_ub,RD_db,T_s] = BL_stall_onset(t,t_ub,t_db,tv0,RD_ub_i,RD_db_i,RD_tv0,theta_i,alpha_lag_i,R_i,upstroke_i,S_i,T_flag_downstroke_i,theta_max_i,theta_min_i,RD_m_i,R_max_i,alpha_lag,q_bar,R,RD,RD_theta,alpha_ds0,alpha_ss,gamma_LS,Tf0);

%% Unsteady breakpoint of separation angles
[alpha1_n,dalpha1_n,alpha1_m,dalpha1_m,alpha1_c,dalpha1_c,dalpha1_db] = BL_alpha_brk(qR,qRdot,R,RD,RD_theta,theta,theta_max,R_max,S,P,T,T_s,T_ub,T_db,upstroke,downstroke_beginning,dalpha1_n_i,dalpha1_db_i,alpha1_0,alpha_ds0,alpha_ss,delta_alpha_0,delta_alpha_1,delta_alpha_2,d_cm,d_cc,z_cm,z_cc);
  
%% Separation points
[f,fprime_n,fprime_m,fprime_c,Sigma2,fprime_n_db,fprime_m_db,fprime_c_db] = BL_sep_points(alpha_bar,alpha_lag,alpha1_n,alpha1_m,alpha1_c,upstroke,downstroke_beginning,theta,theta_min,R,RD,Rdot,S,P,T,fprime_n_i,fprime_m_i,fprime_c_i,fprime_n_db_i,fprime_m_db_i,fprime_c_db_i,T_db,T_s,RD_tv0,T_flag_downstroke,alpha1_0,S1,S2,f0,fb,df0_c,fSig1n,fSig1c,fSig2n,fS2n_ulpr,fS2n_dlpr,fS1n_u,fS1m_u,fS1c_u,fS1n_d,fS1m_d,fS1c_d,fS1n_ud,fS1m_ud,fS1c_ud,fS2n_u,fS2m_u,fS2c_u,fS2n_d,fS2m_d,fS2c_d);

%% Find time of stall onset
[tv0,~,g_v_tv0,f_diff_tv0,RD_tv0,TvL_tv0,theta_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F] = BL_stall_time(tv0,f_diff_tv0,RD_tv0,TvL_tv0,theta_tv0,alpha_lag,alpha_cr,alpha_lag_i,alpha_cr_i,t,t_i,f2prime_n,f,TvL,R,RD,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F,g_v,theta,gamma_TvLs,gamma_TvLc,q_bar);

%% Time delay constants variation
[Tf_n,Tf_m,Tf_c,Ta_theta] = BL_time_constants(in_stall,upstroke,theta,R,RD,P,Ta,Tf0,fSS1,fSS2);

%% State-space
xdot = BLT_state_space(x,xdot,alpha_bar,alpha_tqc,R,A,Tf_n,Tf_m,Tf_c,Ta,fprime_n,fprime_m,fprime_c,Ta_theta);

end