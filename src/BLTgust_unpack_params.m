function [U,b,beta,a_0,a_1,k,A,Ag,Bg,A1,A2,b1,b2,Cg,Cg2,lambda_g,AG,alpha_0L,alpha_ds0,alpha_ss,alpha1_0,gamma_LS,gamma_TvL,delta_alpha_0,delta_alpha_1,delta_alpha_2,kappa,nu_1,nu_2,c_d0,c_m0,c_n_alpha,d_cc,d_cm,df0_c,E0,E1,f0,fb,fSig1n,fSig1c,fSig2n,fS2n_ulpr,fS2n_dlpr,fS1n_u,fS1m_u,fS1c_u,fS1n_d,fS1m_d,fS1c_d,fS1n_ud,fS1m_ud,fS1c_ud,fS2n_u,fS2m_u,fS2c_u,fS2n_d,fS2m_d,fS2c_d,fSS1,fSS2,g_v,K0,K1,K1_f,K2,K3,r0,S1,S2,Ta,Tf0,TvL,Vm,Vn1,Vn2,z_cc,z_cm,gust_profile,gust_options] = BLTgust_unpack_params(params)

%% Test condition
U = params.U;
b = params.b;
beta = params.beta;
a_0 = params.a_0;
a_1 = params.a_1;
k = params.k;
lambda_g = params.lambda_g;

%% State space
A = params.A;
Ag = params.Ag;
Bg = params.Bg;
Cg = params.Cg;
Cg2 = params.Cg2;

%% Airfoil and Mach-dependent
alpha_0L = params.alpha_0L;
alpha_ds0 = params.alpha_ds0;
alpha_ss = params.alpha_ss;
alpha1_0 = params.alpha1_0;
gamma_LS = params.gamma_LS;
gamma_TvL = params.gamma_TvL;
delta_alpha_0 = params.delta_alpha_0; 
delta_alpha_1 = params.delta_alpha_1;
delta_alpha_2 = params.delta_alpha_2;
kappa = params.kappa; 
nu_1 = params.nu_1; 
nu_2 = params.nu_2; 
c_d0 = params.c_d0;  
c_m0 = params.c_m0; 
c_n_alpha = params.c_n_alpha;
d_cc = params.d_cc;
d_cm = params.d_cm;
df0_c = params.df0_c;
E0 = params.E0; 
E1 = params.E1;
f0 = params.f0; 
fb = params.fb;
fSig1n = params.fSig1n;
fSig1c = params.fSig1c;
fSig2n = params.fSig2n;
fS2n_ulpr = params.fS2n_ulpr;
fS2n_dlpr = params.fS2n_dlpr;
fS1n_u = params.fS1n_u;
fS1m_u = params.fS1m_u;
fS1c_u = params.fS1c_u;
fS1n_d = params.fS1n_d;
fS1m_d = params.fS1m_d;
fS1c_d = params.fS1c_d;
fS1n_ud = params.fS1n_ud;
fS1m_ud = params.fS1m_ud;
fS1c_ud = params.fS1c_ud;
fS2n_u = params.fS2n_u;
fS2m_u = params.fS2m_u;
fS2c_u = params.fS2c_u;
fS2n_d = params.fS2n_d;
fS2m_d = params.fS2m_d;
fS2c_d = params.fS2c_d;
fSS1 = params.fSS1; 
fSS2 = params.fSS2; 
g_v = params.g_v;
K0 = params.K0; 
K1 = params.K1; 
K1_f = params.K1_f;
K2 = params.K2; 
K3 = params.K3; 
r0 = params.r0; 
S1 = params.S1; 
S2 = params.S2; 
Ta = params.Ta;
Tf0 = params.Tf0;
TvL = params.TvL;
Vm = params.Vm;
Vn1 = params.Vn1; 
Vn2 = params.Vn2;
z_cc = params.z_cc; 
z_cm = params.z_cm; 

%% Indicial
AG = params.AG;
A1 = params.A1;
A2 = params.A2;
b1 = params.b1;
b2 = params.b2;

%% Gust profile and options
gust_profile = params.gust_profile;
gust_options = params.gust_options;

end