function [c_n,c_nC,c_nf,c_nI,alpha_E,K_f] = BLT_cn_coeff(x,U,b,ah,beta,alphadot,alphaddot,hddot,alpha_tqc,c1,c2,c3,c4,A1,A2,b1,b2,c_n_alpha,f2prime_n,alpha_0L,cn_v)

% % Circulatory thin-airfoil theory: downwash at 3/4-chord and circulatory angle of attack
% w = U*sin(alpha)+b*alphadot;
% alpha_E = 1/U*(1/2*w + c2*c4*(c1+c3)*U^2/b*x(1) + U*(c1*c2+c3*c4)*x(2));

% Effective angle of attack
alpha_E = beta^2*U/b*(A1*b1*x(1)+A2*b2*x(2))+alpha_tqc*(1-A1-A2);

% Circulatory unsteady - attached flow
c_nC = c_n_alpha*sin(alpha_E-alpha_0L);

% Circulatory unsteady - separated flow
K_f = ((1+f2prime_n^(1/2))/2)^2;
c_nf = c_nC*K_f;

% Inertial 
c_nI = pi*b/U^2*(U*alphadot+hddot-ah*b*alphaddot);

% Total: circulatory + inertial + vortex
c_n = c_nf+c_nI+cn_v;

end