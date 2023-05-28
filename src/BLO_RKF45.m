function [outputs,tp,xp,yp,xdotp] = BLO_RKF45(tspan,x0,y0,params,RKF_options)
  
%% Setup the algorithm
% Initialize RK5 algorithm's variables
[delta_b,b_max_it,boundary,display_progress,c5,c45,dt,dt_lim,i,tc,ti,tf,RKF_it_max,RKFtol,tp,xp,yp,xdotp,x_i,xdot_i] = setup_RKF45(tspan,x0,y0,RKF_options);
% Unpack model's parameters
[U,M,b,beta,a_0,a_1,k,A,B,A1,A2,A3,A4,b1,b2,b3,b4,b5,K_a,K_aM,K_q,K_qM,T_I,alpha_0L,alpha1_0,delta_alpha1,eta,kappa,c_d0,c_m0,c_n1,c_n_alpha,Df,E0,f0,fb,K0,K1,K2,S1,S2,Tf0,Tp,Tv0,TvL,x_ac] = BLO_unpack_params(params);
% Initialize model's complementary variables
tv0 = 0; 

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
        [k1,tv0,alpha_i,q_i,alpha1_i] = BLO_dxdt(tc,x_i,xdot_i,tc,tv0,x_i(9),U,b,M,beta,a_0,a_1,k,A,B,A1,A2,b1,b2,K_a,K_q,T_I,alpha1_0,delta_alpha1,c_n1,c_n_alpha,f0,fb,S1,S2,Tf0,Tp,Tv0,TvL);
        x1 = x_i + k1*dt/4;
        [k2,tv0] = BLO_dxdt(tc+dt/4,x1,xdot_i,tc,tv0,x_i(9),U,b,M,beta,a_0,a_1,k,A,B,A1,A2,b1,b2,K_a,K_q,T_I,alpha1_0,delta_alpha1,c_n1,c_n_alpha,f0,fb,S1,S2,Tf0,Tp,Tv0,TvL);
        x2 = x_i + (3/32*k1+9/32*k2)*dt;
        [k3,tv0] = BLO_dxdt(tc+3/8*dt,x2,xdot_i,tc,tv0,x_i(9),U,b,M,beta,a_0,a_1,k,A,B,A1,A2,b1,b2,K_a,K_q,T_I,alpha1_0,delta_alpha1,c_n1,c_n_alpha,f0,fb,S1,S2,Tf0,Tp,Tv0,TvL);
        x3 = x_i + (1932*k1-7200*k2+7296*k3)/2197*dt;
        [k4,tv0] = BLO_dxdt(tc+12/13*dt,x3,xdot_i,tc,tv0,x_i(9),U,b,M,beta,a_0,a_1,k,A,B,A1,A2,b1,b2,K_a,K_q,T_I,alpha1_0,delta_alpha1,c_n1,c_n_alpha,f0,fb,S1,S2,Tf0,Tp,Tv0,TvL);
        x4 = x_i + (439/216*k1-8*k2+3680/513*k3-845/4104*k4)*dt;
        [k5,tv0,alpha_ip1,q_ip1,alpha1_ip1] = BLO_dxdt(tc+dt,x4,xdot_i,tc,tv0,x_i(9),U,b,M,beta,a_0,a_1,k,A,B,A1,A2,b1,b2,K_a,K_q,T_I,alpha1_0,delta_alpha1,c_n1,c_n_alpha,f0,fb,S1,S2,Tf0,Tp,Tv0,TvL);
        x5 = x_i + (-8/27*k1+2*k2-3544/2565*k3+1859/4104*k4-11/40*k5)*dt;
        [k6,tv0] = BLO_dxdt(tc+dt/2,x5,xdot_i,tc,tv0,x_i(9),U,b,M,beta,a_0,a_1,k,A,B,A1,A2,b1,b2,K_a,K_q,T_I,alpha1_0,delta_alpha1,c_n1,c_n_alpha,f0,fb,S1,S2,Tf0,Tp,Tv0,TvL);
        k_step = [k1, k2, k3, k4, k5, k6];
        % Update states rates
        xdot_i = k_step*c5;
        % Update states array
        x_ip1 = x_i + xdot_i*dt;
        % Error between 5th and 4th order RKF approximations   
        eps = norm(k_step*c45)*dt;  
        % Boundary crossing
        boundary = BLO_boundaries(boundary,tc+dt,tc,tv0,alpha_i,alpha1_i,x_i(10),q_i,x_i(9),alpha1_ip1,alpha_ip1,x_ip1(10),q_ip1,x_ip1(9),alpha1_0,c_n1,c_n_alpha,fb,TvL);
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
    % Update storage arrays
    tp(i+1) = tc+dt;      
    xp(:,i+1) = x_ip1;
    xdotp(:,i) = xdot_i;
    yp(:,i+1) = BLO_outputs(tc+dt,x_ip1,tv0,U,b,beta,k,a_0,a_1,M,A1,A2,A3,A4,b1,b2,b3,b4,b5,K_a,K_aM,K_q,K_qM,T_I,alpha_0L,alpha1_0,delta_alpha1,eta,kappa,c_d0,c_m0,c_n1,c_n_alpha,Df,E0,f0,fb,K0,K1,K2,S1,S2,Tf0,Tv0,TvL,x_ac);
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
outputs = BLO_output_vars(tp,xp,yp,xdotp);

end