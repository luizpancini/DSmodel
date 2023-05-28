function [xdot,tv0,alpha,q,alpha1,alpha2,f2prime_tv0] = BLG_dxdt(t,x,xdot,t_i,tv0,so_i,U,b,M,beta,a_0,a_1,k,A,B,A1,A2,b1,b2,K_a,K_q,T_I,alpha1_0,alpha2_0,delta_alpha1,delta_alpha2,c_n1,c_n_alpha,f01,f02,f03,fb1,fb2,fb3,S1,S2,S3,Tf0,Tp,Tv0,TvL,f2prime_tv0)

%% Kinematics
alpha = a_0+a_1*sin(k*U/b*t);
alphadot = a_1*k*U/b*cos(k*U/b*t);
q = 2*alphadot*b/U;
alpha_tqc = alpha+q/2;

%% Rename nonlinear states
c_n_prime = x(9);  
f2prime = x(10);

%% Stall onset criterion
theta = abs(c_n_prime)/c_n1;

%% Unsteady breakpoint of separation angle
[alpha1,alpha2] =  BLG_alpha_brk(alpha,q,alpha1_0,alpha2_0,delta_alpha1,delta_alpha2,f2prime);

%% Separation points
[f,fprime] = BLG_sep_points(c_n_prime,alpha,alpha1_0,alpha1,alpha2_0,alpha2,c_n_alpha,f01,f02,f03,fb1,fb2,fb3,S1,S2,S3);

%% Find time of stall onset
[tv0,tau_v,f2prime_tv0] = BLG_stall_time(tv0,c_n_prime,c_n1,so_i,c_n1,t,t_i,f2prime_tv0,f2prime,f);

%% Time constants modifications
[Tf,Tv] = BLG_time_constants(alpha,q,theta,tau_v,TvL,Tf0,Tv0,alpha2);

%% Normal coefficient variables
[c_nC,c_nP,c_nCdot] = BLO_cn_vars(x,alpha,q,U,b,M,beta,c_n_alpha,A1,A2,b1,b2,K_a,K_q,T_I);

%% Vortex-induced rate of change of vorticity
c_vdot = BLO_vortex_rate(alpha,theta,tau_v,TvL,fprime,f2prime,Tf,c_nC,c_nCdot);

%% State-space
xdot = BLG_state_space(x,xdot,alpha,alpha_tqc,q,A,B,fprime,Tp,Tf,Tv,c_nP,c_vdot);

end