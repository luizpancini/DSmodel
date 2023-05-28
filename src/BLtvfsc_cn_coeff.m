function [c_n,c_nC,c_nf,c_nI,alpha_E,K_f] = BLtvfsc_cn_coeff(U,w_tqc_eff,int_conv_Ia,int_conv_Iq,int_conv_IM,c_n_alpha,f2prime_n,alpha_0L,cn_v)

% Effective angle of attack
alpha_E = w_tqc_eff/U;

% Circulatory unsteady - attached flow
c_nC = c_n_alpha*sin(alpha_E-alpha_0L);

% Circulatory unsteady - separated flow
K_f = ((1+f2prime_n^(1/2))/2)^2;
c_nf = c_nC*K_f;

% Impulsive
c_nI = int_conv_Ia+int_conv_Iq+int_conv_IM;

% Total: circulatory + impulsive + vortex
c_n = c_nf+c_nI+cn_v;

end