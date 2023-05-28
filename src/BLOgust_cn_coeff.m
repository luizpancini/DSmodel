function [c_n,c_nC,c_nf,c_nI,c_nP,c_nv,alpha_Et] = BLOgust_cn_coeff(x,gust_states,U,b,M,beta,Cg,A1,A2,b1,b2,K_a,K_q,T_I,alpha,q,c_n_alpha,f2prime,alpha_0L)

% Effective angles of attack: motion-induced, gust-induced and total
alpha_E = beta^2*U/b*(A1*b1*x(1)+A2*b2*x(2));
alpha_gE = Cg*gust_states;
alpha_Et = alpha_E + alpha_gE;

% Circulatory unsteady - attached flow (gust-augmented)
c_nC = c_n_alpha*(alpha_Et-alpha_0L);

% Impulsive
c_nI = -4/(M*K_a*T_I)*x(3)-1/(M*K_q*T_I)*x(4)+4/M*alpha+1/M*q;

% Total from potential flow
c_nP = c_nC+c_nI;

% Circulatory unsteady - separated flow
K_f = ((1+sqrt(f2prime))/2)^2;
c_nf = c_nC*K_f;

% Vortex
c_nv = x(11);

% Total
c_n = c_nf+c_nI+c_nv;

end

