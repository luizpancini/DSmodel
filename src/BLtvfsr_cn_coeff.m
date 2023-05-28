function [c_n,c_nC,c_nf,c_nI,alpha_E,K_f] = BLtvfsr_cn_coeff(X,U,cnaw_tqc_eff,c_n_alpha,f2prime_n,alpha_0L,cn_v)

% Effective angle of attack
alpha_E = cnaw_tqc_eff/c_n_alpha/U;

% Circulatory unsteady - attached flow
c_nC = c_n_alpha*sin(alpha_E-alpha_0L);

% Circulatory unsteady - separated flow
K_f = ((1+f2prime_n^(1/2))/2)^2;
c_nf = c_nC*K_f;

% Impulsive
c_nI = sum(X(3:5));

% Total: circulatory + impulsive + vortex
c_n = c_nf+c_nI+cn_v;

end
