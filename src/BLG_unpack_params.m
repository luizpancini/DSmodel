function [U,M,b,beta,a_0,a_1,k,A,B,A1,A2,A3,A4,b1,b2,b3,b4,b5,K_a,K_aM,K_q,K_qM,T_I,alpha_0L,alpha1_0,alpha2_0,delta_alpha1,delta_alpha2,eta,kappa,c_d0,c_m0,c_n1,c_n_alpha,Df,E0,Ef,f01,f02,f03,fb1,fb2,fb3,K0,K1,K2,S1,S2,S3,Tf0,Tp,Tv0,TvL,x_ac] = BLG_unpack_params(params)

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

%% Airfoil and Mach-dependent
alpha_0L = params.alpha_0L;
alpha1_0 = params.alpha1_0;
alpha2_0 = params.alpha2_0;
delta_alpha1 = params.delta_alpha1;
delta_alpha2 = params.delta_alpha2;
eta = params.eta;
kappa = params.kappa; 
c_d0 = params.c_d0;  
c_m0 = params.c_m0; 
c_n1 = params.c_n1;
c_n_alpha = params.c_n_alpha;
Df = params.Df;
E0 = params.E0; 
Ef = params.Ef;
f01 = params.f01; 
f02 = params.f02; 
f03 = params.f03; 
fb1 = params.fb1; 
fb2 = params.fb2;
fb3 = params.fb3;
K0 = params.K0; 
K1 = params.K1; 
K2 = params.K2;  
S1 = params.S1; 
S2 = params.S2; 
S3 = params.S3;
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

end