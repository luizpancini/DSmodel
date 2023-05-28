function [theta_min,theta_max,RD_m,R_max,S,upstroke,T_flag_downstroke] = BL_so_ratio_extremes(alpha_lag_i,theta_i,R_i,S_i,alpha_lag_ip1,theta_ip1,R_ip1,q_bar_ip1,RD_ip1,theta_min,theta_max,RD_m,R_max,S,upstroke,T_flag_downstroke)

% Update maximum capped reduced pitch rate
if R_ip1 > R_i
    R_max = R_ip1;
end

% Update maximum stall onset ratio and corresponding delayed capped reduced
% pitch rate (update only for increments in the lagged AoA)
if abs(theta_ip1) > abs(theta_i) && abs(alpha_lag_ip1) > abs(alpha_lag_i)
    theta_max = theta_ip1;
    RD_m = RD_ip1;
end

% Update minimum stall onset ratio (update only for decrements in the lagged AoA)
if abs(theta_ip1) < abs(theta_i) && abs(alpha_lag_ip1) < abs(alpha_lag_i) && theta_ip1*q_bar_ip1>0
    theta_min = theta_ip1;
end

%
S = abs(theta_max) > 1;
upstroke = theta_ip1*q_bar_ip1>=0;
stall_beginning = S & ~S_i;                         
if stall_beginning && ~upstroke
    T_flag_downstroke = false;   
elseif upstroke
    T_flag_downstroke = true;
end

end