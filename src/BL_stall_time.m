function [tv0,tau_v,f_diff_tv0,qR_tv0,R_tv0,RD_tv0,Tv_tv0,theta_tv0,upstroke_tv0,f_diff_tv0_2,RD_tv0_2,upstroke_tv0_2] = BL_stall_time(tv0,f_diff_tv0,qR_tv0,R_tv0,RD_tv0,Tv_tv0,theta_tv0,upstroke_tv0,f_diff_tv0_2,RD_tv0_2,upstroke_tv0_2,t_i,theta_i,t,theta,upstroke,qR,R,RD,f_n,f2prime_n,g_v,Tv)

%% Primary vortex
if abs(theta) >=1 && abs(theta_i) < 1 && t > t_i
    % Linear interpolation for time of vortex shedding 
    tv0 = lininterp1([theta_i theta],[t_i t],1); 
    % Variables at time of vortex shedding
    theta_tv0 = theta;
    f_diff_tv0 = f2prime_n-f_n;
    qR_tv0 = qR;
    R_tv0 = R;
    RD_tv0 = RD;
    upstroke_tv0 = upstroke;  
    Tv_tv0 = Tv;
    % Reset flag to allow secondary vortices
    f_diff_tv0_2 = -1; 
end
% Time since primary vortex shedding 
tau_v = max([0, t-tv0]); 

%% Secondary vortices
if tau_v >= (2+g_v)*Tv && f_diff_tv0_2 == -1
    % Set variables at the time of second vortex
    f_diff_tv0_2 = f2prime_n-f_n;
    RD_tv0_2 = RD;
    upstroke_tv0_2 = upstroke;
end

end