function [c_n,c_nC,c_nf,c_nI,alpha_Et,K_f] = BLgust_cn_coeff(x,M,U,b,ah,beta,c_n_alpha,A1,A2,b1,b2,K_a,K_q,T_I,alpha,alpha_g,alpha_tqc,q,f2prime_n,alpha_0L,cn_v,Cg,AG,gust_states,G,g1,g2)

% Effective angles of attack: pitch-induced, gust-induced and composition
alpha_E = beta^2*U/b*(A1*b1*x(1)+A2*b2*x(2))+alpha_tqc*(1-A1-A2);
alpha_gE = Cg*gust_states+alpha_g*(1-sum(AG));
alpha_Et = alpha_E+alpha_gE;

% Circulatory, gust-augmented unsteady - attached flow
c_nC = c_n_alpha*sin(alpha_Et-alpha_0L)+U/b*beta^2*G*(g1*x(end-2)-g2*x(end-1));

% Circulatory, gust-augmented unsteady - separated flow
K_f = ((1+f2prime_n^(1/2))/2)^2;
c_nf = c_nC*K_f;

% Inertial
c_nI = -4/(M*K_a*T_I)*x(3) - (-2*ah)/(M*K_q*T_I)*x(4) + 4/M*alpha + (-2*ah)/M*q;

% Total: circulatory + inertial + vortex
c_n = c_nf+c_nI+cn_v;

end