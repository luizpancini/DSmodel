function [U,M,b,beta,a_0,a_1,k,A,B,Ag,Bg,Cg,A1,A2,A3,A4,b1,b2,b3,b4,b5,K_a,K_aM,K_q,K_qM,T_I,alpha_0L,alpha1_0,delta_alpha1,eta,kappa,c_d0,c_m0,c_n1,c_n_alpha,Df,E0,f0,fb,K0,K1,K2,S1,S2,Tf0,Tp,Tv0,TvL,x_ac,gust_profile,gust_options] = BLOgust_unpack_params(params)

%% Test condition
U = params.U;
M = params.M;
b = params.b;
a_0 = params.a_0;
a_1 = params.a_1;
k = params.k;
beta = params.beta;

%% State space
A = params.A;
B = params.B;
Ag = params.Ag;
Bg = params.Bg;
Cg = params.Cg;

%% Airfoil and Mach-dependent
alpha_0L = params.alpha_0L;
alpha1_0 = params.alpha1_0;
delta_alpha1 = params.delta_alpha1;
eta = params.eta;
kappa = params.kappa; 
c_d0 = params.c_d0;  
c_m0 = params.c_m0; 
c_n1 = params.c_n1;
c_n_alpha = params.c_n_alpha;
Df = params.Df;
E0 = params.E0; 
f0 = params.f0; 
fb = params.fb; 
K0 = params.K0; 
K1 = params.K1; 
K2 = params.K2;  
S1 = params.S1; 
S2 = params.S2; 
Tf0 = params.Tf0;
Tv0 = params.Tv0;
Tp = params.Tp;
TvL = params.TvL;
x_ac = params.x_ac;

%% Indicial
A1 = params.A1;
A2 = params.A2;
A3 = params.A3;
A4 = params.A4;
b1 = params.b1;
b2 = params.b2;
b3 = params.b3;
b4 = params.b4;
b5 = params.b5;
K_a = params.K_a;
K_aM = params.K_aM;
K_q = params.K_q;
K_qM = params.K_qM;
T_I = params.T_I;

%% Gust
gust_profile = params.gust_profile;
gust_options = params.gust_options;

end