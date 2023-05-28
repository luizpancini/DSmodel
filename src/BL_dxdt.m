function xdot = BL_dxdt(t,x,xdot,dt,tv0,t_ub,t_db,RD_tv0,alpha_lag_i,theta_i,r_i,qR_i,upstroke_i,S_i,P_i,T_flag_downstroke_i,theta_max_i,theta_min_i,RD_m_i,qR_max_i,dalpha1_n_i,dalpha1_db_i,fprime_n_i,fprime_m_i,fprime_c_i,fprime_n_db_i,fprime_m_db_i,fprime_c_db_i,RD_ub_i,RD_db_i,A,U,b,ah,a_0,a_1,k,a_1h,k_h,alpha1_0n,alpha1_0m,alpha1_0c,alpha_ds0,alpha_ss,beta_Sig1n,beta_Sig1c,beta_Sig2n,beta_S2n_lpr,beta_S2c_lpr,beta_S1n_u,beta_S1m_u,beta_S1c_u,beta_S1n_d,beta_S1m_d,beta_S1c_d,beta_S2n_u,beta_S2m_u,beta_S2c_u,beta_S2n_d,beta_S2m_d,beta_S2c_d,gamma_LS,delta_alpha_0,delta_alpha_1,lambda_1,lambda_2,xi,zeta_a,d_cm,d_cc,f0_n,f0_m,f0_c,fb_n,fb_m,fb_c,r0,S1_n,S1_m,S1_c,S2_n,S2_m,S2_c,Ta,Tf,z_cm,z_ccd,z_ccu)

%% Kinematics
[alpha,alpha_plunge,w_QS,Uc,alpha_QS,hdot,w_tqc,q,q_QS,qR,R,r,rdot] = BL_kinematics(t,U,b,ah,k,a_0,a_1,a_1h,k_h,r0,r_i,dt);

%% States
alpha_lag = atan(x(9)/Uc);
RD = x(13);
RD_theta = x(14);

%% Stall onset criterion, stall and motion qualifiers
[alpha_cr,theta,theta_max,theta_min,del_RD_acc,RD_m,qR_max,R_max,S,P,T,upstroke,downstroke_beginning,in_stall,T_flag_downstroke,t_ub,t_db,T_ub,T_db,RD_ub,RD_db,T_s] = BL_stall_onset(t,t_ub,t_db,tv0,RD_ub_i,RD_db_i,RD_tv0,theta_i,alpha_lag_i,qR_i,upstroke_i,S_i,P_i,T_flag_downstroke_i,theta_max_i,theta_min_i,RD_m_i,qR_max_i,alpha_lag,q_QS,qR,R,RD,RD_theta,rdot,alpha_ds0,alpha_ss,gamma_LS,Tf,zeta_a);

%% Unsteady breakpoint of separation angles
[alpha1_n,dalpha1_n,alpha1_m,dalpha1_m,alpha1_c,dalpha1_c] = BL_alpha_brk(qR,R,RD,RD_theta,del_RD_acc,RD_m,theta,theta_min,theta_max,R_max,S,P,T,T_s,T_ub,T_db,upstroke,downstroke_beginning,dalpha1_n_i,dalpha1_db_i,alpha1_0n,alpha1_0m,alpha1_0c,alpha_ds0,alpha_ss,delta_alpha_0,delta_alpha_1,d_cm,d_cc,z_cm,z_ccd,z_ccu);
 
%% Separation points
[f_n,f_m,f_c,fprime_n,fprime_m,fprime_c,Sigma2,fprime_n_db,fprime_m_db,fprime_c_db] = BL_sep_points(alpha_QS,alpha_lag,alpha1_n,alpha1_m,alpha1_c,upstroke,downstroke_beginning,theta,theta_min,R,RD,S,P,T,fprime_n_i,fprime_m_i,fprime_c_i,fprime_n_db_i,fprime_m_db_i,fprime_c_db_i,T_db,T_s,RD_tv0,T_flag_downstroke,alpha1_0n,alpha1_0m,alpha1_0c,beta_Sig1n,beta_Sig1c,beta_Sig2n,beta_S2n_lpr,beta_S2c_lpr,beta_S1n_u,beta_S1m_u,beta_S1c_u,beta_S1n_d,beta_S1m_d,beta_S1c_d,beta_S2n_u,beta_S2m_u,beta_S2c_u,beta_S2n_d,beta_S2m_d,beta_S2c_d,xi,f0_n,f0_m,f0_c,fb_n,fb_m,fb_c,S1_n,S1_m,S1_c,S2_n,S2_m,S2_c);

%% Time delay constants variation
[Tf_n,Tf_m,Tf_c,Ta_theta] = BL_time_constants(in_stall,upstroke,theta,qR,R,RD,P,qR_max,lambda_1,lambda_2,Ta,Tf);

%% State-space
xdot = BL_state_space(x,xdot,w_QS,alpha_QS,w_tqc,q,R,A,fprime_n,fprime_m,fprime_c,Tf_n,Tf_m,Tf_c,Ta,Ta_theta);

end