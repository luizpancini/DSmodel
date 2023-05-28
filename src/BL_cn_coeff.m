function [c_n,c_nC,c_nf,c_nI,alpha_C] = BL_cn_coeff(x,M,U,b,ah,beta,Uc,alpha_bar,q,alpha,w_tqc,hdot,f2prime_n,alpha_0L,c_n_alpha,A1,A2,b1,b2,K_a,K_q,T_I,c_nv)

% Effective downwash
w_E = w_tqc*(1-A1-A2)+beta^2*U/b*(A1*b1*x(1)+A2*b2*x(2));

% Effective circulatory AoA
alpha_C = atan(w_E/Uc);

% Circulatory unsteady - attached flow
c_nC = c_n_alpha/U*(w_E-U*sin(alpha_0L));

% Circulatory unsteady - separated flow
c_nf = c_nC*((1+sqrt(f2prime_n))/2)^2;

% Inertial (at pitch axis)
c_nI = -4/(M*K_a*T_I)*x(3) - (-2*ah)/(M*K_q*T_I)*x(4) + 4/M*alpha_bar + (-2*ah)/M*q;

% Total: circulatory + inertial + vortex
c_n = c_nf+c_nI+c_nv;

end