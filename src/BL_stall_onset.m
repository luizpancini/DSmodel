function [alpha_cr,theta,theta_max,theta_min,del_RD_acc,RD_m,qR_max,R_max,S,P,T,upstroke,downstroke_beginning,in_stall,T_flag_downstroke,t_ub,t_db,T_ub,T_db,RD_ub,RD_db,T_s] = BL_stall_onset(t,t_ub,t_db,tv0,RD_ub,RD_db,RD_tv0,theta_i,alpha_lag_i,qR_i,upstroke_i,S_i,P_i,T_flag_downstroke,theta_max,theta_min,RD_m,qR_max,alpha_lag,q_bar,qR,R,RD,RD_theta,rdot,alpha_ds0,alpha_ss,gamma_LS,Tf,zeta_a)

%% Critical angle for dynamic stall
del_RD_acc = min([zeta_a,abs(rdot)*R*(alpha_lag*q_bar>0)]);
del_RD_acc = 0
alpha_cr = alpha_ss+(alpha_ds0-alpha_ss)*(RD_theta-del_RD_acc);   

%% Dynamic stall onset ratio
theta = alpha_lag/alpha_cr;

%% Update extremes
% Update maximum reduced pitch rate ratio
if qR > qR_i
    qR_max = qR;
end
R_max = min([1,qR_max]);
% Update maximum stall onset ratio and corresponding delayed capped reduced pitch rate (update only for increments in the lagged AoA)
if abs(theta) > abs(theta_i) && abs(alpha_lag) > abs(alpha_lag_i)
    theta_max = theta;
end
% Update minimum stall onset ratio (update only for decrements in the lagged AoA and when in the upstroke)
% if abs(theta) < abs(theta_i) && abs(alpha_lag) < abs(alpha_lag_i) %&& theta*q_bar>0
%     theta_min = theta;
% end

%% Stall condition and motion indicators
% Flags
upstroke = theta*q_bar>=0;                          % TF for upstroke
S = abs(theta_max) > 1;                             % TF for stalled condition (S = 1: stall has occurred, S = 0: stall has not occured)
if S
    P = exp(-gamma_LS*(abs(theta_max)^4-1));          % Indicates how light stall is (P close to 1: light stall, P close to 0: deep stall)
else
    P = P_i;
end
T = exp(-50*RD);                                    % Indicates static condition (T close to 1: static, T close to 0: unsteady)
stall_beginning = S & ~S_i;                         % TF for beginning of stall
in_stall = abs(theta) > 1;                          % TF for currently stalled flow;
if stall_beginning && ~upstroke
    T_flag_downstroke = false;                      % Flag to use T_db instead of T_s for the calculation of separation point in the stalled downstroke
end
if upstroke_i == -1                                 % Neither upstroke nor downstroke (zero pitch rate)
    upstroke_beginning = upstroke;                  % TF for beginning of upstroke
    downstroke_beginning = ~upstroke;               % TF for beginning of downstroke
else
    upstroke_beginning = upstroke & ~upstroke_i;    
    downstroke_beginning = ~upstroke & upstroke_i;  
end
% Time instances
if upstroke_beginning
    t_ub = t;                                       % Time upstroke began
    RD_ub = max([0.1,RD]);                          % Corresponding delayed capped reduced pitch rate
    T_flag_downstroke = true;
    if ~upstroke_i, theta_min = theta; end
end
if downstroke_beginning
    t_db = t;                                       % Time downstroke began
    RD_db = max([0.1,RD]);                          % Corresponding delayed capped reduced pitch rate
    RD_m = RD_theta;
end
% Time decay functions
T_ub = exp(-(t-t_ub)*RD_ub^2/Tf);                  % For build-up of breakpoint angle offset in the upstroke
T_db = exp(-10*(t-t_db)*RD_db^2/Tf);               % For build-up of breakpoint angle offset in the downstroke 
T_s =  exp(-5*(t-tv0)*max([0.1,RD_tv0])^2/Tf);     % For build-up of breakpoint angle offset in the downstroke after stall 

end