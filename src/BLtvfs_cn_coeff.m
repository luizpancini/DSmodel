function [c_n,c_nC,c_nf,c_nI,alpha_E,K_f] = BLtvfs_cn_coeff(x,U,w_tqc,c_n_alpha,f2prime_n,alpha_0L,cn_v)

% Effective downwash at 3/4-chord and angle of attack
w_tqc_eff = w_tqc-sum(x(7:8))/c_n_alpha;

% Effective angle of attack
alpha_E = w_tqc_eff/U;

% Circulatory unsteady - attached flow
c_nC = c_n_alpha*sin(alpha_E-alpha_0L);

% Circulatory unsteady - separated flow
K_f = ((1+f2prime_n^(1/2))/2)^2;
c_nf = c_nC*K_f;

% Inertial
c_nI = sum(x(9:11));

% Total: circulatory + inertial + vortex
c_n = c_nf+c_nI+cn_v;

end