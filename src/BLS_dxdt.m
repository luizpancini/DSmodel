function [xdot,tv0,so_state,so_lim,alpha1,alpha,alphadot,r] = BLS_dxdt(t,x,xdot,t_i,tv0,so_i,so_lim_i,A,B,U,b,a_0,a_1,k,alpha1_0,S1,S2,TvL,Ta,Tf0,r0,alpha_ds0,alpha_ss,f0,fb)

%% Kinematics
[alpha,alphadot,~,alpha_tqc,q,~,r,~,~,~,R,~] = BLO_kinematics(t,U,b,k,a_0,a_1,r0);

%% Stall onset criterion
[so_state,so_lim,theta] = BLS_stall_onset(x(7),R,alpha_ss,alpha_ds0);        

%% alpha1 variation
alpha1 = BLS_alpha_brk(alpha1_0,alpha_ss,alpha_ds0,R,alpha,q);

%% Separation points
[f,f_prime] = BLS_sep_points(so_state,f0,fb,alpha,alpha1,alpha1_0,S1,S2);

%% Find time of stall onset
[tv0,tau_v] = BLS_stall_time(tv0,so_state,so_lim,so_i,so_lim_i,t,t_i);

%% Tf variation
Tf = BLS_time_constants(Tf0,theta,tau_v,TvL,alpha,q,x(8),fb);

%% State-space
xdot = BLS_state_space(x,xdot,alpha,alpha_tqc,q,A,B,f_prime,Tf,Ta,f,Tf0);

end