function [outputs,tp,xp,yp,xdotp] = BLgust_RKF45(tspan,x0,y0,params,RKF_options)

%% Setup the algorithm
% Initialize RK5 algorithm's variables
[delta_b,b_max_it,boundary,display_progress,c5,c45,dt,dt_lim,i,tc,ti,tf,RKF_it_max,RKFtol,tp,xp,yp,xdotp,x_i,xdot_i] = setup_RKF45(tspan,x0,y0,RKF_options);
% Unpack model's parameters
[U,M,b,ah,a_0,a_1,k,beta,A,AG,Ag,Bg,Cg,G,g1,g2,alpha_0L,alpha_ds0,alpha_ss,alpha1_0,gamma_LS,gamma_TvL,delta_alpha_0,delta_alpha_1,delta_alpha_2,kappa,nu_1,nu_2,c_d0,c_m0,c_n_alpha,d_cc,d_cm,df0_c,E0,E1,f0,fb,fSig1n,fSig1c,fSig2n,fS2n_ulpr,fS2n_dlpr,fS1n_u,fS1m_u,fS1c_u,fS1n_d,fS1m_d,fS1c_d,fS1n_ud,fS1m_ud,fS1c_ud,fS2n_u,fS2m_u,fS2c_u,fS2n_d,fS2m_d,fS2c_d,fSS1,fSS2,g_v,K0,K1,K1_f,K2,K3,r0,S1,S2,Ta,Tf0,TvL,Vm,Vn1,Vn2,x_ac,z_cc,z_cm,A1,A2,A3,A4,b1,b2,b3,b4,b5,K_a,K_aM,K_q,K_qM,T_I,gust_profile,gust_options] = BLgust_unpack_params(params);
% Initialize model's complementary variables
tv0 = 0; t_ub = 0; t_db = 0; alpha_lag_i = 0; alpha_cr_i = 1; theta_i = 0; R_i = 0; RD_tv0 = 0; f_diff_tv0 = 0; TvL_tv0 = -1e4; theta_tv0 = 0; g_v_tv0 = g_v; RD_tv0_2 = 1; f_diff_tv0_2 = 0; TvL_tv0_2 = -1e4; V2F = 0; theta_min_i = 0; theta_max_i = 1; RD_m_i = 1; R_max_i = 1;
upstroke_i = -1; dalpha1_n_i = 0; dalpha1_db_i = 0; fprime_n_i = 1; fprime_m_i = 1; fprime_c_i = 1; fprime_n_db_i = 1; fprime_m_db_i = 1; fprime_c_db_i = 1; RD_ub_i = 0; RD_db_i = 0; S_i = 0; T_flag_downstroke_i = true;
% Get gust velocity as a function of time
[wg_fun,wgdot_fun,wgddot_fun] = gust_functions(gust_profile,gust_options,U,b,tf);
% Initialize outputs
yp(:,i) = BLgust_outputs(tc,x_i,tv0,t_ub,t_db,theta_i,alpha_lag_i,R_i,upstroke_i,S_i,T_flag_downstroke_i,dalpha1_n_i,dalpha1_db_i,fprime_n_i,fprime_m_i,fprime_c_i,fprime_n_db_i,fprime_m_db_i,fprime_c_db_i,RD_ub_i,RD_db_i,theta_max_i,theta_min_i,RD_m_i,R_max_i,U,b,ah,beta,k,a_0,a_1,M,x_ac,K0,K1,K1_f,K2,K3,kappa,K_a,K_aM,K_q,K_qM,T_I,b1,b2,b3,b4,b5,A1,A2,A3,A4,AG,Cg,G,g1,g2,alpha_0L,E0,E1,c_m0,c_n_alpha,r0,alpha_ds0,alpha_ss,alpha1_0,S1,S2,Ta,Tf0,f0,fb,Vn1,Vn2,Vm,c_d0,delta_alpha_0,delta_alpha_1,delta_alpha_2,d_cm,d_cc,z_cm,z_cc,nu_1,nu_2,gamma_LS,df0_c,fSig1n,fSig1c,fSig2n,fS2n_ulpr,fS2n_dlpr,fS1n_u,fS1m_u,fS1c_u,fS1n_d,fS1m_d,fS1c_d,fS1n_ud,fS1m_ud,fS1c_ud,fS2n_u,fS2m_u,fS2c_u,fS2n_d,fS2m_d,fS2c_d,fSS1,fSS2,g_v_tv0,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,TvL,wg_fun,wgdot_fun,wgddot_fun);

%% Solve the ODEs            
% Loop over time steps
while tc < tf
    % Reset iterations counters and error
    eps = 10*RKFtol;          
    b_it = 0;                 
    RKF_it = 0;               
    % Loop until convergence is reached
    while (eps > RKFtol || any(boundary < delta_b)) 
        % Adjust for last time step
        dt = min([dt,tf-tc]); 
        % RKF45 steps
        [xdot1,tv0,t_ub,t_db,alpha_bar_i,q_bar_i,alpha_lag_i,alpha_cr_i,alpha1n_i,alpha1m_i,alpha1c_i,dalpha1_n_i,dalpha1_db_i,fprime_n_i,fprime_m_i,fprime_c_i,fprime_n_db_i,fprime_m_db_i,fprime_c_db_i,RD_ub_i,RD_db_i,R_i] = BLgust_dxdt(tc,x_i,xdot_i,tc,tv0,t_ub,t_db,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F,alpha_lag_i,alpha_cr_i,theta_i,R_i,upstroke_i,S_i,T_flag_downstroke_i,theta_max_i,theta_min_i,RD_m_i,R_max_i,dalpha1_n_i,dalpha1_db_i,fprime_n_i,fprime_m_i,fprime_c_i,fprime_n_db_i,fprime_m_db_i,fprime_c_db_i,RD_ub_i,RD_db_i,A,Ag,Bg,U,b,ah,beta,a_0,a_1,k,g1,g2,S1,S2,TvL,Ta,Tf0,r0,alpha_ds0,alpha_ss,alpha1_0,f0,fb,delta_alpha_0,delta_alpha_1,delta_alpha_2,d_cm,d_cc,z_cm,z_cc,gamma_LS,g_v,df0_c,fSS1,fSS2,fSig1n,fSig1c,fSig2n,fS2n_ulpr,fS2n_dlpr,fS1n_u,fS1m_u,fS1c_u,fS1n_d,fS1m_d,fS1c_d,fS1n_ud,fS1m_ud,fS1c_ud,fS2n_u,fS2m_u,fS2c_u,fS2n_d,fS2m_d,fS2c_d,gamma_TvL,wg_fun,wgdot_fun,wgddot_fun); 
        x1 = x_i + xdot1*dt/4;
        [xdot2,tv0,t_ub,t_db] = BLgust_dxdt(tc+dt/4,x1,xdot_i,tc,tv0,t_ub,t_db,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F,alpha_lag_i,alpha_cr_i,theta_i,R_i,upstroke_i,S_i,T_flag_downstroke_i,theta_max_i,theta_min_i,RD_m_i,R_max_i,dalpha1_n_i,dalpha1_db_i,fprime_n_i,fprime_m_i,fprime_c_i,fprime_n_db_i,fprime_m_db_i,fprime_c_db_i,RD_ub_i,RD_db_i,A,Ag,Bg,U,b,ah,beta,a_0,a_1,k,g1,g2,S1,S2,TvL,Ta,Tf0,r0,alpha_ds0,alpha_ss,alpha1_0,f0,fb,delta_alpha_0,delta_alpha_1,delta_alpha_2,d_cm,d_cc,z_cm,z_cc,gamma_LS,g_v,df0_c,fSS1,fSS2,fSig1n,fSig1c,fSig2n,fS2n_ulpr,fS2n_dlpr,fS1n_u,fS1m_u,fS1c_u,fS1n_d,fS1m_d,fS1c_d,fS1n_ud,fS1m_ud,fS1c_ud,fS2n_u,fS2m_u,fS2c_u,fS2n_d,fS2m_d,fS2c_d,gamma_TvL,wg_fun,wgdot_fun,wgddot_fun);
        x2 = x_i + (3/32*xdot1+9/32*xdot2)*dt;
        [xdot3,tv0,t_ub,t_db] = BLgust_dxdt(tc+3/8*dt,x2,xdot_i,tc,tv0,t_ub,t_db,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F,alpha_lag_i,alpha_cr_i,theta_i,R_i,upstroke_i,S_i,T_flag_downstroke_i,theta_max_i,theta_min_i,RD_m_i,R_max_i,dalpha1_n_i,dalpha1_db_i,fprime_n_i,fprime_m_i,fprime_c_i,fprime_n_db_i,fprime_m_db_i,fprime_c_db_i,RD_ub_i,RD_db_i,A,Ag,Bg,U,b,ah,beta,a_0,a_1,k,g1,g2,S1,S2,TvL,Ta,Tf0,r0,alpha_ds0,alpha_ss,alpha1_0,f0,fb,delta_alpha_0,delta_alpha_1,delta_alpha_2,d_cm,d_cc,z_cm,z_cc,gamma_LS,g_v,df0_c,fSS1,fSS2,fSig1n,fSig1c,fSig2n,fS2n_ulpr,fS2n_dlpr,fS1n_u,fS1m_u,fS1c_u,fS1n_d,fS1m_d,fS1c_d,fS1n_ud,fS1m_ud,fS1c_ud,fS2n_u,fS2m_u,fS2c_u,fS2n_d,fS2m_d,fS2c_d,gamma_TvL,wg_fun,wgdot_fun,wgddot_fun);
        x3 = x_i + (1932*xdot1-7200*xdot2+7296*xdot3)/2197*dt;
        [xdot4,tv0,t_ub,t_db] = BLgust_dxdt(tc+12/13*dt,x3,xdot_i,tc,tv0,t_ub,t_db,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F,alpha_lag_i,alpha_cr_i,theta_i,R_i,upstroke_i,S_i,T_flag_downstroke_i,theta_max_i,theta_min_i,RD_m_i,R_max_i,dalpha1_n_i,dalpha1_db_i,fprime_n_i,fprime_m_i,fprime_c_i,fprime_n_db_i,fprime_m_db_i,fprime_c_db_i,RD_ub_i,RD_db_i,A,Ag,Bg,U,b,ah,beta,a_0,a_1,k,g1,g2,S1,S2,TvL,Ta,Tf0,r0,alpha_ds0,alpha_ss,alpha1_0,f0,fb,delta_alpha_0,delta_alpha_1,delta_alpha_2,d_cm,d_cc,z_cm,z_cc,gamma_LS,g_v,df0_c,fSS1,fSS2,fSig1n,fSig1c,fSig2n,fS2n_ulpr,fS2n_dlpr,fS1n_u,fS1m_u,fS1c_u,fS1n_d,fS1m_d,fS1c_d,fS1n_ud,fS1m_ud,fS1c_ud,fS2n_u,fS2m_u,fS2c_u,fS2n_d,fS2m_d,fS2c_d,gamma_TvL,wg_fun,wgdot_fun,wgddot_fun);
        x4 = x_i + (439/216*xdot1-8*xdot2+3680/513*xdot3-845/4104*xdot4)*dt;
        [xdot5,tv0,t_ub,t_db,alpha_bar_ip1,q_bar_ip1,alpha_lag_ip1,alpha_cr_ip1,alpha1n_ip1,alpha1m_ip1,alpha1c_ip1,dalpha1_n_ip1,dalpha1_db_ip1,fprime_n_ip1,fprime_m_ip1,fprime_c_ip1,fprime_n_db_ip1,fprime_m_db_ip1,fprime_c_db_ip1,RD_ub_ip1,RD_db_ip1,R_ip1,g_v_tv0,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F] = BLgust_dxdt(tc+dt,x4,xdot_i,tc,tv0,t_ub,t_db,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F,alpha_lag_i,alpha_cr_i,theta_i,R_i,upstroke_i,S_i,T_flag_downstroke_i,theta_max_i,theta_min_i,RD_m_i,R_max_i,dalpha1_n_i,dalpha1_db_i,fprime_n_i,fprime_m_i,fprime_c_i,fprime_n_db_i,fprime_m_db_i,fprime_c_db_i,RD_ub_i,RD_db_i,A,Ag,Bg,U,b,ah,beta,a_0,a_1,k,g1,g2,S1,S2,TvL,Ta,Tf0,r0,alpha_ds0,alpha_ss,alpha1_0,f0,fb,delta_alpha_0,delta_alpha_1,delta_alpha_2,d_cm,d_cc,z_cm,z_cc,gamma_LS,g_v,df0_c,fSS1,fSS2,fSig1n,fSig1c,fSig2n,fS2n_ulpr,fS2n_dlpr,fS1n_u,fS1m_u,fS1c_u,fS1n_d,fS1m_d,fS1c_d,fS1n_ud,fS1m_ud,fS1c_ud,fS2n_u,fS2m_u,fS2c_u,fS2n_d,fS2m_d,fS2c_d,gamma_TvL,wg_fun,wgdot_fun,wgddot_fun);
        x5 = x_i + (-8/27*xdot1+2*xdot2-3544/2565*xdot3+1859/4104*xdot4-11/40*xdot5)*dt;
        [xdot6,tv0,t_ub,t_db] = BLgust_dxdt(tc+dt/2,x5,xdot_i,tc,tv0,t_ub,t_db,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F,alpha_lag_i,alpha_cr_i,theta_i,R_i,upstroke_i,S_i,T_flag_downstroke_i,theta_max_i,theta_min_i,RD_m_i,R_max_i,dalpha1_n_i,dalpha1_db_i,fprime_n_i,fprime_m_i,fprime_c_i,fprime_n_db_i,fprime_m_db_i,fprime_c_db_i,RD_ub_i,RD_db_i,A,Ag,Bg,U,b,ah,beta,a_0,a_1,k,g1,g2,S1,S2,TvL,Ta,Tf0,r0,alpha_ds0,alpha_ss,alpha1_0,f0,fb,delta_alpha_0,delta_alpha_1,delta_alpha_2,d_cm,d_cc,z_cm,z_cc,gamma_LS,g_v,df0_c,fSS1,fSS2,fSig1n,fSig1c,fSig2n,fS2n_ulpr,fS2n_dlpr,fS1n_u,fS1m_u,fS1c_u,fS1n_d,fS1m_d,fS1c_d,fS1n_ud,fS1m_ud,fS1c_ud,fS2n_u,fS2m_u,fS2c_u,fS2n_d,fS2m_d,fS2c_d,gamma_TvL,wg_fun,wgdot_fun,wgddot_fun);
        xdot_step = [xdot1, xdot2, xdot3, xdot4, xdot5, xdot6];
        % Error between 5th and 4th order RKF approximations   
        eps = norm(xdot_step*c45)*dt;  
        % Boundary crossing
        boundary = BL_boundaries(boundary,tc+dt,tc,tv0,alpha_lag_i,alpha_cr_i,alpha_bar_i,q_bar_i,alpha1n_i,alpha1m_i,alpha1c_i,alpha_lag_ip1,alpha_cr_ip1,alpha_bar_ip1,q_bar_ip1,alpha1n_ip1,alpha1m_ip1,alpha1c_ip1,alpha_ss,TvL,r0);
        % Check boundary crossing tolerance
        if any(boundary < delta_b) 
            % Increase iteration count
            b_it = b_it+1;
            % If maximum number of boundary iterations is reached, break the loop
            if b_it == b_max_it, break; end
            % Halve time step and restart
            dt = dt/2; continue;
        % Check RKF steps tolerance    
        elseif eps > RKFtol
            % Increase iteration count
            RKF_it = RKF_it+1;  
            % If maximum number of RKF iterations is reached, break the loop
            if RKF_it == RKF_it_max, break; end
            % Halve time step and restart
            dt = dt/2; continue;
        end
    end
    % Update states rates 
    xdot_i = xdot_step*c5;
    % Update states array
    x_ip1 = x_i + xdot_i*dt;
    % Update storage arrays
    tp(i+1) = tc+dt;      
    xp(:,i+1) = x_ip1;
    xdotp(:,i) = xdot_i;
    yp(:,i+1) = BLgust_outputs(tc+dt,x_ip1,tv0,t_ub,t_db,theta_i,alpha_lag_i,R_i,upstroke_i,S_i,T_flag_downstroke_i,dalpha1_n_i,dalpha1_db_i,fprime_n_i,fprime_m_i,fprime_c_i,fprime_n_db_i,fprime_m_db_i,fprime_c_db_i,RD_ub_i,RD_db_i,theta_max_i,theta_min_i,RD_m_i,R_max_i,U,b,ah,beta,k,a_0,a_1,M,x_ac,K0,K1,K1_f,K2,K3,kappa,K_a,K_aM,K_q,K_qM,T_I,b1,b2,b3,b4,b5,A1,A2,A3,A4,AG,Cg,G,g1,g2,alpha_0L,E0,E1,c_m0,c_n_alpha,r0,alpha_ds0,alpha_ss,alpha1_0,S1,S2,Ta,Tf0,f0,fb,Vn1,Vn2,Vm,c_d0,delta_alpha_0,delta_alpha_1,delta_alpha_2,d_cm,d_cc,z_cm,z_cc,nu_1,nu_2,gamma_LS,df0_c,fSig1n,fSig1c,fSig2n,fS2n_ulpr,fS2n_dlpr,fS1n_u,fS1m_u,fS1c_u,fS1n_d,fS1m_d,fS1c_d,fS1n_ud,fS1m_ud,fS1c_ud,fS2n_u,fS2m_u,fS2c_u,fS2n_d,fS2m_d,fS2c_d,fSS1,fSS2,g_v_tv0,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,TvL,wg_fun,wgdot_fun,wgddot_fun);
    % Update extremes and flags
    [theta_min_i,theta_max_i,RD_m_i,R_max_i,S_i,upstroke_i,T_flag_downstroke_i] = BL_so_ratio_extremes(x_i(9)+x_i(end),(x_i(9)+x_i(end))/alpha_cr_i,R_i,S_i,x_ip1(9)+x_ip1(end),(x_ip1(9)+x_ip1(end))/alpha_cr_ip1,R_ip1,q_bar_ip1,x_ip1(13),theta_min_i,theta_max_i,RD_m_i,R_max_i,S_i,upstroke_i,T_flag_downstroke_i);
    % Setup next time step
    tc = tc+dt;
    x_i = x_ip1;
    alpha_lag_i = x_ip1(9)+x_ip1(end);
    alpha_cr_i = alpha_cr_ip1;
    R_i = R_ip1;
    theta_i = (x_ip1(9)+x_ip1(end))/alpha_cr_ip1;
    dalpha1_n_i = dalpha1_n_ip1;
    dalpha1_db_i = dalpha1_db_ip1;
    fprime_n_i = fprime_n_ip1;
    fprime_m_i = fprime_m_ip1;
    fprime_c_i = fprime_c_ip1;
    fprime_n_db_i = fprime_n_db_ip1;
    fprime_m_db_i = fprime_m_db_ip1;
    fprime_c_db_i = fprime_c_db_ip1;
    RD_ub_i = RD_ub_ip1;
    RD_db_i = RD_db_ip1;
    i = i+1;
    % Adjust time step
    dt = min(dt_lim,dt*min([4,(eps/RKFtol)^(-1/5)]));
    % Display progress
    if rem(i,1e4) == 0 && display_progress
        disp(['RKF45 progress: ',num2str((tc-ti)/(tf-ti)*100,'%10.2f') '%'])
    end
end

% Assume last time step derivatives equal to previous
xdotp(:,i) = xdotp(:,i-1);

% Truncate pre-allocated arrays
tp = tp(1:i); xp = xp(:,1:i); yp = yp(:,1:i); xdotp = xdotp(:,1:i);

% Get outputs structure
outputs = BLgust_output_vars(tp,xp,yp,xdotp);

end