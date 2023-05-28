function [xdot,tv0,alpha,q,alpha1] = BLOgust_dxdt(t,x,xdot,t_i,tv0,so_i,U,b,M,beta,a_0,a_1,k,A,B,Ag,Bg,Cg,A1,A2,b1,b2,K_a,K_q,T_I,alpha_0L,alpha1_0,delta_alpha1,c_n1,c_n_alpha,f0,fb,S1,S2,Tf0,Tp,Tv0,TvL,wg_fun,wgdot_fun)

%% Kinematics
% Gust-induced AoA and its rate
alpha_g = get_alpha_g(wg_fun,wgdot_fun,t,U);
% Motion-induced AoA and its rate, and nondimensional pitch rate
alpha = a_0 + a_1*sin(k*U/b*t);
alphadot = a_1*k*U/b*cos(k*U/b*t);
q = 2*alphadot*b/U;
% Angle of attack at 3/4-chord
alpha_tqc = alpha+q/2;

%% Nonlinear and gust states
so_state = x(9);  
f2prime = x(10);
gust_states = x(13:end);
gust_states_rates = Ag*gust_states+Bg*alpha_g;

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
Tf = 4*Tf0;

%% Normal coefficient variables
[c_nC,c_nP,c_nCdot] = BLOgust_cn_vars(x,alpha,q,U,b,M,beta,c_n_alpha,alpha_0L,A1,A2,b1,b2,K_a,K_q,T_I,Cg,gust_states,gust_states_rates);

%% Vortex-induced rate of change of vorticity
c_vdot = BLO_vortex_rate(alpha,theta,tau_v,TvL,fprime,f2prime,Tf,c_nC,c_nCdot);

%% State-space
xdot = BLOgust_state_space(x,xdot,alpha,alpha_tqc,q,A,B,f,fprime,Tp,Tf,Tv,Tf0,c_nP,c_vdot,gust_states_rates);

end