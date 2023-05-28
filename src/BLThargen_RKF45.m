function [outputs,tp,xp,yp,xdotp] = BLThargen_RKF45(tspan,x0,y0,params,RKF_options)

%% Setup the algorithm
% Initialize RKF algorithm's variables
[display_progress,dtf_steps,c5,c45,dt,dt_min,dt_max,i,tc,ti,tf,RKF_it_max,RKFtol,tp,xp,yp,xdotp,x_i,xdot_i,xdot] = setup_RKF45(tspan,x0,y0,RKF_options);
% Unpack model's parameters
[airfoil,b,ah,a_inf,a_0,a_1,k,a_1h,k_h,psi_ha,update_params,U_0,U_1,k_U,psi_Ua,dh,d_0,d_1,k_f,psi_fa,gust_profile,lambda_g,gust_options,gust_ind,A1,A2,b1,b2,A1W,A2W,b1W,b2W,AG,bG,Cg,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,alpha_0L,alpha1_0n,alpha1_0m,alpha1_0c,alpha_ds0,alpha_ss,beta_Sig1n,beta_Sig1c,beta_Sig2n,beta_S2n_lpr,beta_S2c_lpr,beta_S1n_u,beta_S1m_u,beta_S1c_u,beta_S1n_d,beta_S1m_d,beta_S1c_d,beta_S2n_u,beta_S2m_u,beta_S2c_u,beta_S2n_d,beta_S2m_d,beta_S2c_d,gamma_LS,delta_alpha_0,delta_alpha_1,eps_fn,eps_fm,eps_fh,eta,kappa_0,kappa_1,kappa_2,kappa_3,lambda_1,lambda_2,mu_v2,nu_1,nu_2,nu_3,nu_4,nu_5,chi_u,chi_d,xi,zeta_a,c_d0,c_m0,c_n_alpha,d_cc,d_cm,E0,E1,f0_n,f0_m,f0_c,fb_n,fb_m,fb_c,g_v,g_v2,K0,K1,K2,r0,S1_n,S1_m,S1_c,S2_n,S2_m,S2_c,Ta,Tf,Tv,Tv2,Vn1,Vn2,Vn3,Vm,Vc,z_ccd,z_ccu,z_cm] = BLThargen_unpack_params(params);
% Initialize model's complementary variables
tv0 = -1e9; t_ub = 0; t_db = 0; alpha_lag_i = 0; theta_i = 0; r_i = 1; qR_i = 1; qR_tv0 = 1; R_tv0 = 1; RD_tv0 = 1; f_diff_tv0 = 0; upstroke_tv0 = 0; Tv_tv0 = -1e4; theta_tv0 = 0; RD_tv0_2 = 1; f_diff_tv0_2 = -1; upstroke_tv0_2 = 0; theta_min_i = 0; theta_max_i = 1; RD_m_i = 1; qR_max_i = 1;
upstroke_i = -1; dalpha1_n_i = 0; dalpha1_db_i = 0; fprime_n_i = 1; fprime_m_i = 1; fprime_c_i = 1; fprime_n_db_i = 1; fprime_m_db_i = 1; fprime_c_db_i = 1; RD_ub_i = 0; RD_db_i = 0; S_i = 0; P_i = 0; T_flag_downstroke_i = true;
c_n_alpha_i = c_n_alpha; F_theta0_i = 0; F2_theta0_i = 0;
% Get gust velocity as a function of time
[wg_fun,wgdot_fun,wgddot_fun] = gust_functions(gust_profile,gust_options,U_0,b,tf);
% Set initial outputs
[yp(:,i),F_theta0_i,F2_theta0_i] = BLThargen_outputs(tc,-inf,x_i,tv0,t_ub,t_db,F_theta0_i,F2_theta0_i,r_i,qR_i,theta_i,alpha_lag_i,upstroke_i,S_i,P_i,T_flag_downstroke_i,dalpha1_n_i,dalpha1_db_i,fprime_n_i,fprime_m_i,fprime_c_i,fprime_n_db_i,fprime_m_db_i,fprime_c_db_i,RD_ub_i,RD_db_i,theta_max_i,theta_min_i,RD_m_i,qR_max_i,c_n_alpha_i,qR_tv0,R_tv0,RD_tv0,f_diff_tv0,theta_tv0,upstroke_tv0,Tv_tv0,f_diff_tv0_2,RD_tv0_2,upstroke_tv0_2,airfoil,b,ah,dh,a_inf,a_0,a_1,k,a_1h,k_h,psi_ha,U_0,U_1,k_U,psi_Ua,d_0,d_1,k_f,psi_fa,lambda_g,AG,bG,Cg,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,alpha_0L,alpha1_0n,alpha1_0m,alpha1_0c,alpha_ds0,alpha_ss,beta_Sig1n,beta_Sig1c,beta_Sig2n,beta_S2n_lpr,beta_S2c_lpr,beta_S1n_u,beta_S1m_u,beta_S1c_u,beta_S1n_d,beta_S1m_d,beta_S1c_d,beta_S2n_u,beta_S2m_u,beta_S2c_u,beta_S2n_d,beta_S2m_d,beta_S2c_d,gamma_LS,delta_alpha_0,delta_alpha_1,eps_fn,eps_fm,eps_fh,eta,kappa_0,kappa_1,kappa_2,kappa_3,lambda_1,lambda_2,mu_v2,nu_1,nu_2,nu_3,nu_4,nu_5,chi_u,chi_d,xi,zeta_a,c_d0,c_m0,c_n_alpha,d_cc,d_cm,E0,E1,f0_n,f0_m,f0_c,fb_n,fb_m,fb_c,g_v,g_v2,K0,K1,K2,r0,S1_n,S1_m,S1_c,S2_n,S2_m,S2_c,Ta,Tf,Tv,Tv2,Vn1,Vn2,Vn3,Vm,Vc,z_ccd,z_ccu,z_cm,wg_fun,wgdot_fun,wgddot_fun);

%% Solve the ODEs            
% Loop over time steps
while tc < tf
    % Reset iterations counters and error
    eps = 10*RKFtol;                          
    RKF_it = 0;   
    if tc >= 0.0924*3
        here=1;
    end
    % Loop until convergence is reached
    while eps > RKFtol 
        % Adjust for last time step
        dt = min([dt,tf-tc]); 
        % RKF45 steps        
        for RKF_step=1:6
            switch RKF_step
                case 1, x_step = x_i;
                case 2, x_step = x_i + (xdot(:,1)/4)*dt;
                case 3, x_step = x_i + (3/32*xdot(:,1)+9/32*xdot(:,2))*dt;
                case 4, x_step = x_i + (1932*xdot(:,1)-7200*xdot(:,2)+7296*xdot(:,3))/2197*dt;
                case 5, x_step = x_i + (439/216*xdot(:,1)-8*xdot(:,2)+3680/513*xdot(:,3)-845/4104*xdot(:,4))*dt;
                case 6, x_step = x_i + (-8/27*xdot(:,1)+2*xdot(:,2)-3544/2565*xdot(:,3)+1859/4104*xdot(:,4)-11/40*xdot(:,5))*dt;    
            end
            t_step = tc + dtf_steps(RKF_step)*dt;
            xdot(:,RKF_step) = BLThargen_dxdt(t_step,x_step,xdot_i,tc,tv0,t_ub,t_db,RD_tv0,alpha_lag_i,theta_i,r_i,qR_i,upstroke_i,S_i,P_i,T_flag_downstroke_i,theta_max_i,theta_min_i,RD_m_i,qR_max_i,dalpha1_n_i,dalpha1_db_i,fprime_n_i,fprime_m_i,fprime_c_i,fprime_n_db_i,fprime_m_db_i,fprime_c_db_i,RD_ub_i,RD_db_i,c_n_alpha_i,airfoil,a_inf,b,ah,a_0,a_1,k,a_1h,k_h,psi_ha,U_0,U_1,k_U,psi_Ua,d_0,d_1,k_f,psi_fa,lambda_g,A1,A2,b1,b2,A1W,A2W,b1W,b2W,bG,T10,T11,alpha1_0n,alpha1_0m,alpha1_0c,alpha_ds0,alpha_ss,beta_Sig1n,beta_Sig1c,beta_Sig2n,beta_S2n_lpr,beta_S2c_lpr,beta_S1n_u,beta_S1m_u,beta_S1c_u,beta_S1n_d,beta_S1m_d,beta_S1c_d,beta_S2n_u,beta_S2m_u,beta_S2c_u,beta_S2n_d,beta_S2m_d,beta_S2c_d,gamma_LS,delta_alpha_0,delta_alpha_1,lambda_1,lambda_2,xi,zeta_a,d_cm,d_cc,f0_n,f0_m,f0_c,fb_n,fb_m,fb_c,r0,S1_n,S1_m,S1_c,S2_n,S2_m,S2_c,Ta,Tf,z_cm,z_ccd,z_ccu,wg_fun);
        end
        % Error between 5th and 4th order RKF approximations   
        eps = norm(xdot*c45)*dt;  
        % Check RKF steps tolerance    
        if eps > RKFtol 
            % Increase iteration count
            RKF_it = RKF_it+1;  
            % If maximum number of RKF iterations is reached, break the loop
            if RKF_it == RKF_it_max || dt/2 < dt_min, break; end
            % Halve time step and restart
            dt = dt/2; continue;
        end
    end
    % Update states rates 
    xdot_i = xdot*c5;
    % Update states array
    x_ip1 = x_i+xdot_i*dt;
    % Update storage arrays
    tp(i+1) = tc+dt;      
    xp(:,i+1) = x_ip1;
    xdotp(:,i) = xdot_i;
    [yp(:,i+1),F_theta0_i,F2_theta0_i,U,t_ub,t_db,r_i,qR_i,theta_i,alpha_lag_i,dalpha1_n_i,fprime_n_i,fprime_m_i,fprime_c_i,fprime_n_db_i,fprime_m_db_i,fprime_c_db_i,RD_ub_i,RD_db_i,theta_min_i,theta_max_i,RD_m_i,qR_max_i,S_i,P_i,upstroke_i,T_flag_downstroke_i,c_n_alpha_i,tv0,f_diff_tv0,qR_tv0,R_tv0,RD_tv0,Tv_tv0,theta_tv0,upstroke_tv0,f_diff_tv0_2,RD_tv0_2,upstroke_tv0_2] = ...
    BLThargen_outputs(tc+dt,tc,x_ip1,tv0,t_ub,t_db,F_theta0_i,F2_theta0_i,r_i,qR_i,theta_i,alpha_lag_i,upstroke_i,S_i,P_i,T_flag_downstroke_i,dalpha1_n_i,dalpha1_db_i,fprime_n_i,fprime_m_i,fprime_c_i,fprime_n_db_i,fprime_m_db_i,fprime_c_db_i,RD_ub_i,RD_db_i,theta_max_i,theta_min_i,RD_m_i,qR_max_i,c_n_alpha_i,qR_tv0,R_tv0,RD_tv0,f_diff_tv0,theta_tv0,upstroke_tv0,Tv_tv0,f_diff_tv0_2,RD_tv0_2,upstroke_tv0_2,airfoil,b,ah,dh,a_inf,a_0,a_1,k,a_1h,k_h,psi_ha,U_0,U_1,k_U,psi_Ua,d_0,d_1,k_f,psi_fa,lambda_g,AG,bG,Cg,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,alpha_0L,alpha1_0n,alpha1_0m,alpha1_0c,alpha_ds0,alpha_ss,beta_Sig1n,beta_Sig1c,beta_Sig2n,beta_S2n_lpr,beta_S2c_lpr,beta_S1n_u,beta_S1m_u,beta_S1c_u,beta_S1n_d,beta_S1m_d,beta_S1c_d,beta_S2n_u,beta_S2m_u,beta_S2c_u,beta_S2n_d,beta_S2m_d,beta_S2c_d,gamma_LS,delta_alpha_0,delta_alpha_1,eps_fn,eps_fm,eps_fh,eta,kappa_0,kappa_1,kappa_2,kappa_3,lambda_1,lambda_2,mu_v2,nu_1,nu_2,nu_3,nu_4,nu_5,chi_u,chi_d,xi,zeta_a,c_d0,c_m0,c_n_alpha,d_cc,d_cm,E0,E1,f0_n,f0_m,f0_c,fb_n,fb_m,fb_c,g_v,g_v2,K0,K1,K2,r0,S1_n,S1_m,S1_c,S2_n,S2_m,S2_c,Ta,Tf,Tv,Tv2,Vn1,Vn2,Vn3,Vm,Vc,z_ccd,z_ccu,z_cm,wg_fun,wgdot_fun,wgddot_fun);
    % Setup next time step
    tc = tc+dt;
    x_i = x_ip1;
    i = i+1;
    % Update Mach-dependent airfoil parameters
    if any(U_1) && rem(i,update_params) == 0 
        [~,alpha_0L,alpha_ds0,alpha_ss,alpha1_0n,alpha1_0m,alpha1_0c,beta_Sig1n,beta_Sig1c,beta_Sig2n,beta_S2n_lpr,beta_S2c_lpr,beta_S1n_u,beta_S1m_u,beta_S1c_u,beta_S1n_d,beta_S1m_d,beta_S1c_d,beta_S2n_u,beta_S2m_u,beta_S2c_u,beta_S2n_d,beta_S2m_d,beta_S2c_d,gamma_LS,delta_alpha_0,delta_alpha_1,eps_fn,eps_fm,eps_fh,eta,kappa_0,kappa_1,kappa_2,kappa_3,lambda_1,lambda_2,mu_v2,nu_1,nu_2,nu_3,nu_4,nu_5,chi_u,chi_d,xi,zeta_a,c_d0,c_m0,c_n_alpha,d_cc,d_cm,E0,E1,f0_n,f0_m,f0_c,fb_n,fb_m,fb_c,g_v,g_v2,K0,K1,K2,r0,S1_n,S1_m,S1_c,S2_n,S2_m,S2_c,Ta,Tf,Tv,Tv2,Vn1,Vn2,Vn3,Vm,Vc,z_ccd,z_ccu,z_cm] = BL_airfoil_parameters(struct(),airfoil,U/a_inf,U,b);
    end
    % Adjust time step
    dt = max([dt_min,min(dt_max,dt*min([4,(eps/RKFtol)^(-1/5)]))]);
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
outputs = BLhargen_output_vars(tp,xp,yp,xdotp);

end