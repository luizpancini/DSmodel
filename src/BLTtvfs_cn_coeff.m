function [c_n,c_nC,c_nf,c_nI,alpha_E,K_f] = BLTtvfs_cn_coeff(x,b,U,Udot,alpha,alphadot,alphaddot,w_tqc,c_n_alpha,f2prime_n,alpha_0L,cn_v)

% Effective angle of attack
alpha_E = 1/U*(w_tqc-sum(x(7:8)));

% Circulatory unsteady - attached flow
c_nC = c_n_alpha*sin(alpha_E-alpha_0L);

% Circulatory unsteady - separated flow
K_f = ((1+f2prime_n^(1/2))/2)^2;
c_nf = c_nC*K_f;

% Inertial 
c_nI = pi*b/U^2*(U*alphadot+alpha*Udot+1/2*b*alphaddot);

% Total: circulatory + inertial + vortex
c_n = c_nf+c_nI+cn_v;

end