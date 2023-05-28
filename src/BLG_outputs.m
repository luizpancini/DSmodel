function y = BLG_outputs(t,x,tv0,U,b,beta,k,a_0,a_1,M,A1,A2,A3,A4,b1,b2,b3,b4,b5,K_a,K_aM,K_q,K_qM,T_I,alpha_0L,alpha1_0,alpha2_0,delta_alpha1,delta_alpha2,eta,kappa,c_d0,c_m0,c_n1,c_n_alpha,Df,E0,Ef,f01,f02,f03,fb1,fb2,fb3,K0,K1,K2,S1,S2,S3,Tf0,Tv0,TvL,x_ac,f2prime_tv0)

%% AoA and non-dimensional pitch rate
alpha = a_0+a_1*sin(k*U/b*t);
alpha_dot = a_1*k*U/b*cos(k*U/b*t);
q = 2*alpha_dot*b/U;

%% Rename nonlinear states
c_n_prime = x(9);  
f2prime = x(10);

%% Stall onset criterion
theta = abs(c_n_prime)/c_n1;

%% Unsteady breakpoint of separation angle
[alpha1,alpha2,dalpha1,dalpha2] =  BLG_alpha_brk(alpha,q,alpha1_0,alpha2_0,delta_alpha1,delta_alpha2,f2prime);

%% Separation points
[f,fprime,alpha_lag] = BLG_sep_points(c_n_prime,alpha,alpha1_0,alpha1,alpha2_0,alpha2,c_n_alpha,f01,f02,f03,fb1,fb2,fb3,S1,S2,S3);

%% Find time of stall onset
tau_v = max([0 t-tv0]);

%% Time constants modifications
[Tf,Tv] = BLG_time_constants(alpha,q,theta,tau_v,TvL,Tf0,Tv0,alpha2);

%% Airload coefficients
% Normal
[c_n,c_nC,c_nf,c_nI,c_nP,c_nv,alpha_E] = BLO_cn_coeff(x,U,b,M,beta,A1,A2,b1,b2,K_a,K_q,T_I,alpha,q,c_n_alpha,f2prime,alpha_0L);
% Pitching moment
[c_m,c_mC,c_mf,c_mI,c_mv,dCP] = BLG_cm_coeff(x,U,b,M,beta,c_nC,c_nf,c_nv,alpha,q,A3,A4,b3,b4,b5,K_aM,K_qM,T_I,K0,K1,K2,kappa,tau_v,TvL,c_m0,x_ac,alpha_lag,alpha1_0,alpha2_0,alpha2,f,fprime,f2prime,f02,fb2,S2);
% Chordwise
c_c = BLG_cc_coeff(alpha,c_n_prime,f,f2prime,f2prime_tv0,theta,alpha_E,alpha_0L,eta,c_d0,c_n1,c_n_alpha,Df,E0,Ef);

%% Outputs
y = [alpha; alpha1; alpha_E; theta; tau_v; c_c; c_m; c_mC; c_mf; c_mI; c_mv; c_n; c_nC; c_nf; c_nI; c_nP; c_nv; dalpha1; dCP; f; fprime; q; Tf; Tv];

end