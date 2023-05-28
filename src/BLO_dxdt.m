function [xdot,tv0,alpha,q,alpha1] = BLO_dxdt(t,x,xdot,t_i,tv0,so_i,U,b,M,beta,a_0,a_1,k,A,B,A1,A2,b1,b2,K_a,K_q,T_I,alpha1_0,delta_alpha1,c_n1,c_n_alpha,f0,fb,S1,S2,Tf0,Tp,Tv0,TvL)

%% Kinematics
alpha = a_0+a_1*sin(k*U/b*t);
alphadot = a_1*k*U/b*cos(k*U/b*t);
q = 2*alphadot*b/U;
alpha_tqc = alpha+q/2;

%% Rename nonlinear states
so_state = x(9);  
f2prime = x(10);

%% Stall onset criterion
theta = abs(so_state)/c_n1;

%% Unsteady breakpoint of separation angle
alpha1 =  BLO_alpha_brk(alpha,q,alpha1_0,delta_alpha1,f2prime);

%% Separation points
[f,fprime] = BLO_sep_points(so_state,alpha,alpha1_0,alpha1,f0,fb,S1,S2,c_n_alpha);

%% Find time of stall onset
[tv0,tau_v] = BLS_stall_time(tv0,so_state,c_n1,so_i,c_n1,t,t_i);

%% Time constants modifications
[Tf,Tv] = BLO_time_constants(alpha,q,theta,tau_v,TvL,Tf0,Tv0,f2prime,fb);

%% Normal coefficient variables
[c_nC,c_nP,c_nCdot] = BLO_cn_vars(x,alpha,q,U,b,M,beta,c_n_alpha,A1,A2,b1,b2,K_a,K_q,T_I);

%% Vortex-induced rate of change of vorticity
c_vdot = BLO_vortex_rate(alpha,theta,tau_v,TvL,fprime,f2prime,Tf,c_nC,c_nCdot);

%% State-space
xdot = BLO_state_space(x,xdot,alpha,alpha_tqc,q,A,B,f,fprime,Tp,Tf,Tv,Tf0,c_nP,c_vdot);

end