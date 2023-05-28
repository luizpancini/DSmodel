function [U,M,b,ah,a_0,a_1,k,a_1h,k_h,beta,alpha_0L,alpha1_0n,alpha1_0m,alpha1_0c,alpha_ds0,alpha_ss,beta_Sig1n,beta_Sig1c,beta_Sig2n,beta_S2n_lpr,beta_S2c_lpr,beta_S1n_u,beta_S1m_u,beta_S1c_u,beta_S1n_d,beta_S1m_d,beta_S1c_d,beta_S2n_u,beta_S2m_u,beta_S2c_u,beta_S2n_d,beta_S2m_d,beta_S2c_d,gamma_LS,delta_alpha_0,delta_alpha_1,eta,kappa_0,kappa_1,kappa_2,kappa_3,lambda_1,lambda_2,mu_v2,nu_1,nu_2,nu_3,nu_4,nu_5,chi_u,chi_d,xi,zeta_a,c_d0,c_m0,c_n_alpha,d_cc,d_cm,E0,E1,f0_n,f0_m,f0_c,fb_n,fb_m,fb_c,g_v,g_v2,K0,K1,K2,r0,S1_n,S1_m,S1_c,S2_n,S2_m,S2_c,Ta,Tf,Tv,Tv2,Vn1,Vn2,Vn3,Vm,Vc,z_ccd,z_ccu,z_cm,A,A1,A2,A3,A4,b1,b2,b3,b4,b5,K_a,K_aM,K_q,K_qM,T_I] = BL_unpack_params(params)

%% Test condition
U = params.U;
M = params.M;
b = params.b;
ah = params.ah;
a_0 = params.a_0;
a_1 = params.a_1;
k = params.k;
beta = params.beta;
a_1h = params.a_1h;
k_h = params.k_h;

%% Airfoil and Mach-dependent
alpha_0L = params.alpha_0L;
alpha_ds0 = params.alpha_ds0;
alpha_ss = params.alpha_ss;
alpha1_0n = params.alpha1_0n;
alpha1_0m = params.alpha1_0m;
alpha1_0c = params.alpha1_0c;
beta_Sig1n = params.beta_Sig1n;
beta_Sig1c = params.beta_Sig1c;
beta_Sig2n = params.beta_Sig2n;
beta_S2n_lpr = params.beta_S2n_lpr;
beta_S2c_lpr = params.beta_S2c_lpr;
beta_S1n_u = params.beta_S1n_u;
beta_S1m_u = params.beta_S1m_u;
beta_S1c_u = params.beta_S1c_u;
beta_S1n_d = params.beta_S1n_d;
beta_S1m_d = params.beta_S1m_d;
beta_S1c_d = params.beta_S1c_d;
beta_S2n_u = params.beta_S2n_u;
beta_S2m_u = params.beta_S2m_u;
beta_S2c_u = params.beta_S2c_u;
beta_S2n_d = params.beta_S2n_d;
beta_S2m_d = params.beta_S2m_d;
beta_S2c_d = params.beta_S2c_d;
gamma_LS = params.gamma_LS;
delta_alpha_0 = params.delta_alpha_0; 
delta_alpha_1 = params.delta_alpha_1;
eta = params.eta;
kappa_0 = params.kappa_0; 
kappa_1 = params.kappa_1; 
kappa_2 = params.kappa_2; 
kappa_3 = params.kappa_3; 
lambda_1 = params.lambda_1; 
lambda_2 = params.lambda_2; 
mu_v2 = params.mu_v2;
nu_1 = params.nu_1; 
nu_2 = params.nu_2;
nu_3 = params.nu_3;
nu_4 = params.nu_4;
nu_5 = params.nu_5;
chi_u = params.chi_u;
chi_d = params.chi_d;
xi = params.xi;
zeta_a = params.zeta_a;
c_d0 = params.c_d0;  
c_m0 = params.c_m0; 
c_n_alpha = params.c_n_alpha;
d_cc = params.d_cc;
d_cm = params.d_cm;
E0 = params.E0; 
E1 = params.E1;
f0_n = params.f0_n; 
f0_m = params.f0_m;
f0_c = params.f0_c;
fb_n = params.fb_n;
fb_m = params.fb_m;
fb_c = params.fb_c;
g_v = params.g_v;
g_v2 = params.g_v2;
K0 = params.K0; 
K1 = params.K1; 
K2 = params.K2;  
r0 = params.r0; 
S1_n = params.S1_n;
S1_m = params.S1_m; 
S1_c = params.S1_c; 
S2_n = params.S2_n; 
S2_m = params.S2_m; 
S2_c = params.S2_c; 
Ta = params.Ta;
Tf = params.Tf;
Tv = params.Tv;
Tv2 = params.Tv2;
Vm = params.Vm;
Vc = params.Vc;
Vn1 = params.Vn1; 
Vn2 = params.Vn2;
Vn3 = params.Vn3;
z_ccd = params.z_ccd; 
z_ccu = params.z_ccu;
z_cm = params.z_cm; 

%% Indicial
A = params.A;
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