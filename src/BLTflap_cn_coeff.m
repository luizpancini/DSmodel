function [c_n,c_nC,c_nf,c_nI,alpha_Et,K_f] = BLTflap_cn_coeff(x,U,b,beta,alpha,alphadot,alphaddot,alpha_tqc,deltadot,deltaddot,delta_qs,T1,T4,eps_n,c1,c2,c3,c4,A1,A2,b1,b2,c_n_alpha,f2prime_n,alpha_0L,cn_v)

% % Circulatory thin-airfoil theory: downwash at 3/4-chord and circulatory angle of attack
% w = U*sin(alpha)+b*alphadot;
% alpha_E = 1/U*(1/2*w + c2*c4*(c1+c3)*U^2/b*x(1) + U*(c1*c2+c3*c4)*x(2));

% Effective angle of attack: airfoil-induced, flap-induced and total
alpha_E = beta^2*U/b*(A1*b1*x(1)+A2*b2*x(2))+alpha_tqc*(1-A1-A2);
alpha_Ef = beta^2*U/b*(A1*b1*x(9)+A2*b2*x(10))+delta_qs*(1-A1-A2);
alpha_Et = alpha_E+eps_n*alpha_Ef;

% Circulatory unsteady - attached flow
c_nC = c_n_alpha*sin(alpha_Et-alpha_0L);

% Circulatory unsteady - separated flow
K_f = ((1+f2prime_n^(1/2))/2)^2;
c_nf = c_nC*K_f;

% Impulsive: airfoil-induced, flap-induced and total
c_nIa = pi*b/U^2*(U*alphadot+1/2*b*alphaddot);
c_nIf = -eps_n*b/U^2*(U*T4*deltadot+b*T1*deltaddot);
c_nI = c_nIa+c_nIf;

% Total: circulatory + impulsive + vortex
c_n = c_nf+c_nI+cn_v;

end
