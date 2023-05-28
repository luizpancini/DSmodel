function y = BLO_outputs(t,x,tv0,U,b,beta,k,a_0,a_1,M,A1,A2,A3,A4,b1,b2,b3,b4,b5,K_a,K_aM,K_q,K_qM,T_I,alpha_0L,alpha1_0,delta_alpha1,eta,kappa,c_d0,c_m0,c_n1,c_n_alpha,Df,E0,f0,fb,K0,K1,K2,S1,S2,Tf0,Tv0,TvL,x_ac)

%% Kinematics
alpha = a_0+a_1*sin(k*U/b*t);
alphadot = a_1*k*U/b*cos(k*U/b*t);
q = 2*alphadot*b/U;

%% Rename nonlinear states
so_state = x(9);  
f2prime = x(10);

%% Stall onset criterion
theta = abs(so_state)/c_n1;

%% Unsteady breakpoint of separation angle
[alpha1,dalpha1] =  BLO_alpha_brk(alpha,q,alpha1_0,delta_alpha1,f2prime);

%% Separation points
[f,fprime] = BLO_sep_points(so_state,alpha,alpha1_0,alpha1,f0,fb,S1,S2,c_n_alpha);

%% Find time of stall onset
tau_v = max([0 t-tv0]);

%% Time constants modifications
[Tf,Tv] = BLO_time_constants(alpha,q,theta,tau_v,TvL,Tf0,Tv0,f2prime,fb);

%% Airload coefficients
% Normal
[c_n,c_nC,c_nf,c_nI,c_nP,c_nv,alpha_E] = BLO_cn_coeff(x,U,b,M,beta,A1,A2,b1,b2,K_a,K_q,T_I,alpha,q,c_n_alpha,f2prime,alpha_0L);
% Pitching moment
[c_m,c_mC,c_mf,c_mI,c_mv,dCP] = BLO_cm_coeff(x,U,b,M,beta,c_nC,c_nf,c_nv,alpha,q,A3,A4,b3,b4,b5,K_aM,K_qM,T_I,K0,K1,K2,kappa,tau_v,TvL,c_m0,x_ac);
% Chordwise
c_c = BLO_cc_coeff(alpha,so_state,f2prime,theta,alpha_E,alpha_0L,eta,c_d0,c_n1,c_n_alpha,Df,E0);

%% Outputs
y = [alpha; alpha1; alpha_E; theta; tau_v; c_c; c_m; c_mC; c_mf; c_mI; c_mv; c_n; c_nC; c_nf; c_nI; c_nP; c_nv; dalpha1; dCP; f; fprime; q; Tf; Tv];

end