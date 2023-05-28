function [f_n,f_m,f_c,fprime_n,fprime_m,fprime_c,sigma,fprime_n_db,fprime_m_db,fprime_c_db] = BL_sep_points(alpha_QS,alpha_lag,alpha1_n,alpha1_m,alpha1_c,upstroke,downstroke_beginning,theta,theta_min,R,RD,S,P,T,fprime_n_i,fprime_m_i,fprime_c_i,fprime_n_db,fprime_m_db,fprime_c_db,T_db,T_s,RD_tv0,T_flag_downstroke,alpha1_0n,alpha1_0m,alpha1_0c,beta_Sig1n,beta_Sig1c,beta_Sig2n,beta_S2n_lpr,beta_S2c_lpr,beta_S1n_u,beta_S1m_u,beta_S1c_u,beta_S1n_d,beta_S1m_d,beta_S1c_d,beta_S2n_u,beta_S2m_u,beta_S2c_u,beta_S2n_d,beta_S2m_d,beta_S2c_d,xi,f0_n,f0_m,f0_c,fb_n,fb_m,fb_c,S1_n,S1_m,S1_c,S2_n,S2_m,S2_c)

%% Light stall and in-and-out-of-stall factors
psi = P^3*(1-T);                                           % Factor defining the flattening of f_prime for light stall, (1-T) to not apply on static conditions
sigma = min([0.5,abs(theta_min)^2.5])*(theta_min*theta>0); % Factor defining how close to stalling angle has been at beginning of upstroke
xi_RD_tv0 = xi*((1-cos(2*pi*RD_tv0))/2);

%% Decay time functions
if T_flag_downstroke 
    T_decay_n = T_db^3; T_decay_m = T_db^4; T_decay_c = T_db^3;
else
    T_decay_n = T_s^3;  T_decay_m = T_s^4;  T_decay_c = T_s^3;
end

%% Separation points based on lagged AoA at downstroke's or stall in downstroke beginning
if downstroke_beginning || ~T_flag_downstroke
    fprime_n_db = fprime_n_i; 
    fprime_m_db = fprime_m_i; 
    fprime_c_db = fprime_c_i; 
end

%% Normal force unsteady lagged separation point based on lagged AoA - fprime_n
if abs(alpha_lag)<=alpha1_n
    alpha_diff = abs(alpha_lag)-alpha1_n;
    if upstroke  
        S1_prime = S1_n*(1+beta_S1n_u*RD);   
        fprime_n = 1-(1-fb_n)*exp(alpha_diff/S1_prime);
    else                
        S1_prime = S1_n*(1+beta_S1n_d*RD); 
%         fprime_n = fprime_n_db*T_decay_n+(1-(1-fb_n)*exp(alpha_diff/S1_prime))*(1-T_decay_n);
        fprime_n = 1-(1-fb_n)*exp(alpha_diff/S1_prime);
    end
else
    alpha_diff = alpha1_n-abs(alpha_lag); 
    if upstroke  
        S2_prime = S2_n*(1+beta_S2n_u*RD+beta_S2n_lpr*(1-RD)^2*(1-T)+beta_Sig2n*sigma*RD*(1-RD)); % beta_S2n_ulpr*(1-RD^(1/2)) to delay the dettachment at low pitch rates, so as +beta_Sig2n*Sigma2*RD*(1-RD) [RD*(1-RD) to have greater effect on medium pitch rates, e.g. frames 10117 and 10118], -D_S2n to match rapid c_n drop and slow growth at stall in almost static conditions
        fprime_n = f0_n+(fb_n-f0_n)*exp(alpha_diff/S2_prime);
    else                
        % If did not stall, do not allow f' to enter in stalled conditions
        % in the downstroke by using the extended formula for alpha < alpha_brk - 
        % check results on e.g. frames 8203 and 9106 
        if ~S, f0_n = fb_n-0.2; end 
%         if ~S           % Unstalled
%             S1_prime = beta_S1n_ud*S1_n; % Flatten the reattachment
%             fprime_n = max([f0_n,fprime_n_db*T_decay_n+(1-(1-fb_n)*exp(-alpha_diff/S1_prime))*(1-T_decay_n)]);
%         else            % Stalled
            f0_n = min([0.25,(1-sqrt(R))*(f0_n+xi_RD_tv0*sigma)]);
            S2_prime = S2_n*(1+beta_S2n_d*RD+beta_S2n_lpr*(1-RD)^2*(1-T)+beta_Sig2n*sigma*RD*(1-RD)+beta_Sig1n*psi*R); 
%             fprime_n = fprime_n_db*T_decay_n+(f0_n+(fb_n-f0_n)*exp(alpha_diff/S2_prime))*(1-T_decay_n);
            fprime_n = f0_n+(fb_n-f0_n)*exp(alpha_diff/S2_prime);
%         end
    end
end

%% Pitching moment unsteady lagged separation point based on lagged AoA - fprime_m
if abs(alpha_lag)<=alpha1_m
    alpha_diff = abs(alpha_lag)-alpha1_m;
    if upstroke  
        S1_prime = S1_m*(1+beta_S1m_u*RD);   
        fprime_m = 1-(1-fb_m)*exp(alpha_diff/S1_prime);
    else                
        S1_prime = S1_m*(1+beta_S1m_d*RD); 
%         fprime_m = fprime_m_db*T_decay_m+(1-(1-fb_m)*exp(alpha_diff/S1_prime))*(1-T_decay_m);
        fprime_m = 1-(1-fb_m)*exp(alpha_diff/S1_prime);
    end   
else
    alpha_diff = alpha1_m-abs(alpha_lag);
    if upstroke  
        S2_prime = S2_m*(1+beta_S2m_u*RD); 
        fprime_m = f0_m+(fb_m-f0_m)*exp(alpha_diff/S2_prime);
    else                
        if ~S, f0_m = fb_m-0.2; end
%         if ~S           % Unstalled - use extended formula for alpha < alpha_brk
%             S1_prime = beta_S1m_ud*S1_m; 
%             fprime_m =  max([f0_m,fprime_m_db*T_decay_m+(1-(1-fb_m)*exp(-alpha_diff/S1_prime))*(1-T_decay_m)]);
%         else            % Stalled
            f0_m = min([0.25,(1-sqrt(R))*(f0_m+xi_RD_tv0*sigma)]);
            S2_prime = S2_m*(1+beta_S2m_d*RD); 
%             fprime_m = fprime_m_db*T_decay_m+(f0_m+(fb_m-f0_m)*exp(alpha_diff/S2_prime))*(1-T_decay_m);
            fprime_m = f0_m+(fb_m-f0_m)*exp(alpha_diff/S2_prime);
%         end
    end
end

%% Chordwise force unsteady lagged separation point based on lagged AoA - fprime_c
if abs(alpha_lag)<=alpha1_c
    alpha_diff = abs(alpha_lag)-alpha1_c;
    if upstroke  
        S1_prime = S1_c*(1+beta_S1c_u*RD); 
        fprime_c = 1-(1-fb_c)*exp(alpha_diff/S1_prime);
    else                
        S1_prime = S1_c*(1+beta_S1c_d*RD);     
%         fprime_c = fprime_c_db*T_decay_c+(1-(1-fb_c)*exp(alpha_diff/S1_prime))*(1-T_decay_c);
        fprime_c = 1-(1-fb_c)*exp(alpha_diff/S1_prime);
    end
else
    alpha_diff = alpha1_c-abs(alpha_lag);
    if upstroke  
        S2_prime = S2_c*(1+beta_S2c_u*RD+beta_S2c_lpr*(1-RD)^2*(1-T)); 
        fprime_c = f0_c+(fb_c-f0_c)*exp(alpha_diff/S2_prime);
    else  
        if ~S, f0_c = fb_c-0.2; end 
%         if ~S           % Unstalled - use extended formula for alpha < alpha_brk
%             S1_prime = beta_S1c_ud*S1_c; 
%             fprime_c =  max([f0_c,fprime_c_db*T_decay_c+(1-(1-fb_c)*exp(-alpha_diff/S1_prime))*(1-T_decay_c)]);
%         else            % Stalled
            f0_c = min([0.25,(1-sqrt(R))*(f0_c+xi_RD_tv0*sigma)]);
            S2_prime = S2_c*(1+beta_S2c_d*RD+beta_S2c_lpr*(1-RD)^2*(1-T)+beta_Sig1c*psi*R); 
%             fprime_c = fprime_c_db*T_decay_c+(f0_c+(fb_c-f0_c)*exp(alpha_diff/S2_prime))*(1-T_decay_c);
            fprime_c = f0_c+(fb_c-f0_c)*exp(alpha_diff/S2_prime);
%         end
    end
end

%% Quasi-steady separation point based on quasi-steady AoA - f
% Normal
if abs(alpha_QS)<=alpha1_0n
    f_n = 1-(1-fb_n)*exp((abs(alpha_QS)-alpha1_0n)/S1_n);
else
    f_n = f0_n+(fb_n-f0_n)*exp((alpha1_0n-abs(alpha_QS))/S2_n);
end
% Pitching moment
if abs(alpha_QS)<=alpha1_0m
    f_m = 1-(1-fb_m)*exp((abs(alpha_QS)-alpha1_0m)/S1_m);
else
    f_m = f0_m+(fb_m-f0_m)*exp((alpha1_0m-abs(alpha_QS))/S2_m);
end
% Chordwise
if abs(alpha_QS)<=alpha1_0c
    f_c = 1-(1-fb_c)*exp((abs(alpha_QS)-alpha1_0c)/S1_c);
else
    f_c = f0_c+(fb_c-f0_c)*exp((alpha1_0c-abs(alpha_QS))/S2_c);
end

end