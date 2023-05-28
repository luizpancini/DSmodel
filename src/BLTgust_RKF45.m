function [outputs,tp,xp,yp,xdotp] = BLTgust_RKF45(tspan,x0,y0,params,RKF_options)
  
%% Setup the algorithm
% Initialize RK5 algorithm's variables
[delta_b,b_max_it,boundary,display_progress,c5,c45,dt,dt_lim,i,tc,ti,tf,RKF_it_max,RKFtol,tp,xp,yp,xdotp,x_i,xdot_i] = setup_RKF45(tspan,x0,y0,RKF_options);
% Unpack model'S parameters
[U,b,beta,a_0,a_1,k,A,Ag,Bg,A1,A2,b1,b2,Cg,Cg2,lambda_g,AG,alpha_0L,alpha_ds0,alpha_ss,alpha1_0,gamma_LS,gamma_TvL,delta_alpha_0,delta_alpha_1,delta_alpha_2,kappa,nu_1,nu_2,c_d0,c_m0,c_n_alpha,d_cc,d_cm,df0_c,E0,E1,f0,fb,fSig1n,fSig1c,fSig2n,fS2n_ulpr,fS2n_dlpr,fS1n_u,fS1m_u,fS1c_u,fS1n_d,fS1m_d,fS1c_d,fS1n_ud,fS1m_ud,fS1c_ud,fS2n_u,fS2m_u,fS2c_u,fS2n_d,fS2m_d,fS2c_d,fSS1,fSS2,g_v,K0,K1,K1_f,K2,K3,r0,S1,S2,Ta,Tf0,TvL,Vm,Vn1,Vn2,z_cc,z_cm,gust_profile,gust_options] = BLTgust_unpack_params(params);
% Initialize model'S complementary variables
tv0 = 0; so_i = 0; so_lim_i = 1; RD_tv0 = 0; f_diff_tv0 = 0; TvL_tv0 = -1e4; theta_tv0 = 0; RD_tv0_2 = 1; f_diff_tv0_2 = 0; TvL_tv0_2 = -1e4; theta_min = 0; theta_max = 1; RD_m = 1; V2F = 0;
% Get gust velocity as a function of time and set initial outputs
[wg_fun,wgdot_fun] = gust_functions(gust_profile,gust_options,U,b,tf);
[yp(:,i),F_theta0,F2_theta0] = BLTgust_outputs(tc,-inf,x_i,nan,nan,tv0,U,b,beta,k,a_0,a_1,A1,A2,b1,b2,Cg,Cg2,lambda_g,AG,K0,K1,K1_f,K2,K3,kappa,alpha_0L,E0,E1,c_m0,c_n_alpha,r0,alpha_ds0,alpha_ss,alpha1_0,S1,S2,Ta,Tf0,f0,fb,Vn1,Vn2,Vm,c_d0,delta_alpha_0,delta_alpha_1,delta_alpha_2,d_cm,d_cc,z_cm,z_cc,nu_1,nu_2,gamma_LS,df0_c,fSig1n,fSig1c,fSig2n,fS2n_ulpr,fS2n_dlpr,fS1n_u,fS1m_u,fS1c_u,fS1n_d,fS1m_d,fS1c_d,fS1n_ud,fS1m_ud,fS1c_ud,fS2n_u,fS2m_u,fS2c_u,fS2n_d,fS2m_d,fS2c_d,fSS1,fSS2,g_v,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,theta_max,theta_min,RD_m,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,TvL,wg_fun,wgdot_fun);

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
        [k1,tv0,so_i,so_lim_i,alpha1n_i,alpha1m_i,alpha1c_i,alpha_bar_i,q_bar_i] = BLTgust_dxdt(tc,x_i,xdot_i,tc,tv0,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F,so_i,so_lim_i,A,Ag,Bg,U,b,a_0,a_1,k,lambda_g,S1,S2,TvL,Ta,Tf0,r0,alpha_ds0,alpha_ss,alpha1_0,f0,fb,delta_alpha_0,delta_alpha_1,delta_alpha_2,d_cm,d_cc,z_cm,z_cc,gamma_LS,g_v,df0_c,fSS1,fSS2,fSig1n,fSig1c,fSig2n,fS2n_ulpr,fS2n_dlpr,fS1n_u,fS1m_u,fS1c_u,fS1n_d,fS1m_d,fS1c_d,fS1n_ud,fS1m_ud,fS1c_ud,fS2n_u,fS2m_u,fS2c_u,fS2n_d,fS2m_d,fS2c_d,gamma_TvL,theta_max,theta_min,RD_m,wg_fun,wgdot_fun); 
        x1 = x_i + k1*dt/4;
        [k2,tv0] = BLTgust_dxdt(tc+dt/4,x1,xdot_i,tc,tv0,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F,so_i,so_lim_i,A,Ag,Bg,U,b,a_0,a_1,k,lambda_g,S1,S2,TvL,Ta,Tf0,r0,alpha_ds0,alpha_ss,alpha1_0,f0,fb,delta_alpha_0,delta_alpha_1,delta_alpha_2,d_cm,d_cc,z_cm,z_cc,gamma_LS,g_v,df0_c,fSS1,fSS2,fSig1n,fSig1c,fSig2n,fS2n_ulpr,fS2n_dlpr,fS1n_u,fS1m_u,fS1c_u,fS1n_d,fS1m_d,fS1c_d,fS1n_ud,fS1m_ud,fS1c_ud,fS2n_u,fS2m_u,fS2c_u,fS2n_d,fS2m_d,fS2c_d,gamma_TvL,theta_max,theta_min,RD_m,wg_fun,wgdot_fun);
        x2 = x_i + (3/32*k1+9/32*k2)*dt;
        [k3,tv0] = BLTgust_dxdt(tc+3/8*dt,x2,xdot_i,tc,tv0,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F,so_i,so_lim_i,A,Ag,Bg,U,b,a_0,a_1,k,lambda_g,S1,S2,TvL,Ta,Tf0,r0,alpha_ds0,alpha_ss,alpha1_0,f0,fb,delta_alpha_0,delta_alpha_1,delta_alpha_2,d_cm,d_cc,z_cm,z_cc,gamma_LS,g_v,df0_c,fSS1,fSS2,fSig1n,fSig1c,fSig2n,fS2n_ulpr,fS2n_dlpr,fS1n_u,fS1m_u,fS1c_u,fS1n_d,fS1m_d,fS1c_d,fS1n_ud,fS1m_ud,fS1c_ud,fS2n_u,fS2m_u,fS2c_u,fS2n_d,fS2m_d,fS2c_d,gamma_TvL,theta_max,theta_min,RD_m,wg_fun,wgdot_fun);
        x3 = x_i + (1932*k1-7200*k2+7296*k3)/2197*dt;
        [k4,tv0] = BLTgust_dxdt(tc+12/13*dt,x3,xdot_i,tc,tv0,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F,so_i,so_lim_i,A,Ag,Bg,U,b,a_0,a_1,k,lambda_g,S1,S2,TvL,Ta,Tf0,r0,alpha_ds0,alpha_ss,alpha1_0,f0,fb,delta_alpha_0,delta_alpha_1,delta_alpha_2,d_cm,d_cc,z_cm,z_cc,gamma_LS,g_v,df0_c,fSS1,fSS2,fSig1n,fSig1c,fSig2n,fS2n_ulpr,fS2n_dlpr,fS1n_u,fS1m_u,fS1c_u,fS1n_d,fS1m_d,fS1c_d,fS1n_ud,fS1m_ud,fS1c_ud,fS2n_u,fS2m_u,fS2c_u,fS2n_d,fS2m_d,fS2c_d,gamma_TvL,theta_max,theta_min,RD_m,wg_fun,wgdot_fun);
        x4 = x_i + (439/216*k1-8*k2+3680/513*k3-845/4104*k4)*dt;
        [k5,tv0,so_ip1,so_lim_ip1,alpha1n_ip1,alpha1m_ip1,alpha1c_ip1,alpha_bar_ip1,q_bar_ip1,g_v_tv0,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F] = BLTgust_dxdt(tc+dt,x4,xdot_i,tc,tv0,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F,so_i,so_lim_i,A,Ag,Bg,U,b,a_0,a_1,k,lambda_g,S1,S2,TvL,Ta,Tf0,r0,alpha_ds0,alpha_ss,alpha1_0,f0,fb,delta_alpha_0,delta_alpha_1,delta_alpha_2,d_cm,d_cc,z_cm,z_cc,gamma_LS,g_v,df0_c,fSS1,fSS2,fSig1n,fSig1c,fSig2n,fS2n_ulpr,fS2n_dlpr,fS1n_u,fS1m_u,fS1c_u,fS1n_d,fS1m_d,fS1c_d,fS1n_ud,fS1m_ud,fS1c_ud,fS2n_u,fS2m_u,fS2c_u,fS2n_d,fS2m_d,fS2c_d,gamma_TvL,theta_max,theta_min,RD_m,wg_fun,wgdot_fun);
        x5 = x_i + (-8/27*k1+2*k2-3544/2565*k3+1859/4104*k4-11/40*k5)*dt;
        [k6,tv0] = BLTgust_dxdt(tc+dt/2,x5,xdot_i,tc,tv0,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,V2F,so_i,so_lim_i,A,Ag,Bg,U,b,a_0,a_1,k,lambda_g,S1,S2,TvL,Ta,Tf0,r0,alpha_ds0,alpha_ss,alpha1_0,f0,fb,delta_alpha_0,delta_alpha_1,delta_alpha_2,d_cm,d_cc,z_cm,z_cc,gamma_LS,g_v,df0_c,fSS1,fSS2,fSig1n,fSig1c,fSig2n,fS2n_ulpr,fS2n_dlpr,fS1n_u,fS1m_u,fS1c_u,fS1n_d,fS1m_d,fS1c_d,fS1n_ud,fS1m_ud,fS1c_ud,fS2n_u,fS2m_u,fS2c_u,fS2n_d,fS2m_d,fS2c_d,gamma_TvL,theta_max,theta_min,RD_m,wg_fun,wgdot_fun);
        k_step = [k1, k2, k3, k4, k5, k6];
        % Update states rates 
        xdot_i = k_step*c5;
        % Update minima and maxima of stall onset ratio
        [theta_min,theta_max,RD_m] = BL_so_ratio_extremes(q_bar_i,so_i,so_ip1,so_lim_i,so_lim_ip1,x_i(7),theta_min,theta_max,RD_m);
        % Error between 5th and 4th order RKF approximations   
        eps = norm(k_step*c45)*dt;  
        % Boundary crossing
        boundary = BL_boundaries(boundary,tc+dt,tc,tv0,so_i,so_lim_i,alpha_bar_i,q_bar_i,alpha1n_i,alpha1m_i,alpha1c_i,so_ip1,so_lim_ip1,alpha_bar_ip1,q_bar_ip1,alpha1n_ip1,alpha1m_ip1,alpha1c_ip1,alpha_ss,TvL,r0);
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
    % Update states array
    x_ip1 = x_i + xdot_i*dt;
    % Update storage arrays
    tp(i+1) = tc+dt;      
    xp(:,i+1) = x_ip1;
    xdotp(:,i) = xdot_i;
    [yp(:,i+1),F_theta0,F2_theta0] = BLTgust_outputs(tc+dt,dt,x_ip1,F_theta0,F2_theta0,tv0,U,b,beta,k,a_0,a_1,A1,A2,b1,b2,Cg,Cg2,lambda_g,AG,K0,K1,K1_f,K2,K3,kappa,alpha_0L,E0,E1,c_m0,c_n_alpha,r0,alpha_ds0,alpha_ss,alpha1_0,S1,S2,Ta,Tf0,f0,fb,Vn1,Vn2,Vm,c_d0,delta_alpha_0,delta_alpha_1,delta_alpha_2,d_cm,d_cc,z_cm,z_cc,nu_1,nu_2,gamma_LS,df0_c,fSig1n,fSig1c,fSig2n,fS2n_ulpr,fS2n_dlpr,fS1n_u,fS1m_u,fS1c_u,fS1n_d,fS1m_d,fS1c_d,fS1n_ud,fS1m_ud,fS1c_ud,fS2n_u,fS2m_u,fS2c_u,fS2n_d,fS2m_d,fS2c_d,fSS1,fSS2,g_v_tv0,RD_tv0,f_diff_tv0,TvL_tv0,theta_tv0,theta_max,theta_min,RD_m,f_diff_tv0_2,RD_tv0_2,TvL_tv0_2,TvL,wg_fun,wgdot_fun);
    % Setup next time step
    tc = tc+dt;
    x_i = x_ip1;
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
outputs = BLTgust_output_vars(tp,xp,yp,xdotp);

end