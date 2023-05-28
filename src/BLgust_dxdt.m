function [xdot,tv0,t_ub,t_db,alpha_bar,q_bar,alpha_lag,alpha_cr,alpha1_n,alpha1_m,alpha1_c,dalpha1_n,dalpha1_db,fprime_n,fprime_m,fprime_c,fprime_n_db,fprime_m_db,fprime_c_db,RD_ub,RD_db,R,g_v_tv0,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F] = ...
          BLgust_dxdt(t,x,xdot,t_i,tv0,t_ub,t_db,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F,alpha_lag_i,alpha_cr_i,theta_i,R_i,upstroke_i,S_i,T_flag_downstroke_i,theta_max_i,theta_min_i,RD_m_i,R_max_i,dalpha1_n_i,dalpha1_db_i,fprime_n_i,fprime_m_i,fprime_c_i,fprime_n_db_i,fprime_m_db_i,fprime_c_db_i,RD_ub_i,RD_db_i,A,Ag,Bg,U,b,ah,beta,a_0,a_1,k,g1,g2,S1,S2,TvL,Ta,Tf0,r0,alpha_ds0,alpha_ss,alpha1_0,f0,fb,delta_alpha_0,delta_alpha_1,delta_alpha_2,d_cm,d_cc,z_cm,z_cc,gamma_LS,g_v,df0_c,fSS1,fSS2,fSig1n,fSig1c,fSig2n,fS2n_ulpr,fS2n_dlpr,fS1n_u,fS1m_u,fS1c_u,fS1n_d,fS1m_d,fS1c_d,fS1n_ud,fS1m_ud,fS1c_ud,fS2n_u,fS2m_u,fS2c_u,fS2n_d,fS2m_d,fS2c_d,gamma_TvL,wg_fun,wgdot_fun,wgddot_fun)

%% Kinematics
[alpha_g,alpha,alphadot,alphaddot,alpha_tqc,q,qdot,r,rdot,qR,qRdot,R,Rdot,alpha_bar,q_bar] = gust_kinematics(t,U,b,ah,k,a_0,a_1,r0,wg_fun,wgdot_fun,wgddot_fun);

%% States
alpha_lag = x(9)+x(end); % Composition of pitch-induced and gust-induced AoA lagged
f2prime_n = x(10);
RD = x(13);
RD_theta = x(14);
gust_states = x(15:end-3);

%% Stall onset criterion, stall and motion qualifiers
[alpha_cr,theta,theta_max,theta_min,RD_m,R_max,S,P,T,upstroke,downstroke_beginning,in_stall,T_flag_downstroke,t_ub,t_db,T_ub,T_db,RD_ub,RD_db,T_s] = BL_stall_onset(t,t_ub,t_db,tv0,RD_ub_i,RD_db_i,RD_tv0,theta_i,alpha_lag_i,R_i,upstroke_i,S_i,T_flag_downstroke_i,theta_max_i,theta_min_i,RD_m_i,R_max_i,alpha_lag,q_bar,R,RD,RD_theta,alpha_ds0,alpha_ss,gamma_LS,Tf0);

%% Unsteady breakpoint of separation angles
[alpha1_n,dalpha1_n,alpha1_m,dalpha1_m,alpha1_c,dalpha1_c,dalpha1_db] = BL_alpha_brk(qR,qRdot,R,RD,RD_theta,theta,theta_max,R_max,S,P,T,T_s,T_ub,T_db,upstroke,downstroke_beginning,dalpha1_n_i,dalpha1_db_i,alpha1_0,alpha_ds0,alpha_ss,delta_alpha_0,delta_alpha_1,delta_alpha_2,d_cm,d_cc,z_cm,z_cc);
  
%% Separation points
[f,fprime_n,fprime_m,fprime_c,Sigma2,fprime_n_db,fprime_m_db,fprime_c_db] = BL_sep_points(alpha_bar,alpha_lag,alpha1_n,alpha1_m,alpha1_c,upstroke,downstroke_beginning,theta,theta_min,R,RD,Rdot,S,P,T,fprime_n_i,fprime_m_i,fprime_c_i,fprime_n_db_i,fprime_m_db_i,fprime_c_db_i,T_db,T_s,RD_tv0,T_flag_downstroke,alpha1_0,S1,S2,f0,fb,df0_c,fSig1n,fSig1c,fSig2n,fS2n_ulpr,fS2n_dlpr,fS1n_u,fS1m_u,fS1c_u,fS1n_d,fS1m_d,fS1c_d,fS1n_ud,fS1m_ud,fS1c_ud,fS2n_u,fS2m_u,fS2c_u,fS2n_d,fS2m_d,fS2c_d);

%% Find time of stall onset
[tv0,~,g_v_tv0,f_diff_tv0,RD_tv0,TvL_tv0,theta_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F] = BL_stall_time(tv0,f_diff_tv0,RD_tv0,TvL_tv0,theta_tv0,alpha_lag,alpha_cr,alpha_lag_i,alpha_cr_i,t,t_i,f2prime_n,f,TvL,R,RD,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F,g_v,theta,gamma_TvL,q_bar);

%% Time delay constants variation
[Tf_n,Tf_m,Tf_c,Ta_theta] = BL_time_constants(in_stall,upstroke,theta,R,RD,Ta,Tf0,fSS1,fSS2);

%% State-space
xdot = BLgust_state_space(x,xdot,alpha,alpha_tqc,q,R,A,Tf_n,Tf_m,Tf_c,Ta,fprime_n,fprime_m,fprime_c,Ta_theta,Ag,Bg,gust_states,alpha_g,U,b,beta,g1,g2);

end