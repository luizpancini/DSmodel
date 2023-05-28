function [c_nC,c_nP,c_nCdot] = BLOgust_cn_vars(x,alpha,q,U,b,M,beta,c_n_alpha,alpha_0L,A1,A2,b1,b2,K_a,K_q,T_I,Cg,gust_states,gust_states_rates)

% Effective angles of attack: motion-induced, gust-induced and total
alpha_E = beta^2*U/b*(A1*b1*x(1)+A2*b2*x(2));
alpha_gE = Cg*gust_states;
alpha_Et = alpha_E + alpha_gE;

% Circulatory unsteady - attached flow (gust-augmented)
c_nC = c_n_alpha*(alpha_Et-alpha_0L);

% Circulatory normal coefficient rate
c_nCdot = c_n_alpha*U/b*beta^2*(A1*b1*(-U/b*beta^2*b1*x(1)+alpha+q/2)+A2*b2*(-U/b*beta^2*b2*x(2)+alpha+q/2)) + c_n_alpha*Cg*gust_states_rates;

% Impulsive normal coefficient
c_nI = - 4/(M*K_a*T_I)*x(3)-1/(M*K_q*T_I)*x(4)+4/M*alpha+1/M*q;

% Total potential flow normal coefficient
c_nP = c_nC + c_nI;

end