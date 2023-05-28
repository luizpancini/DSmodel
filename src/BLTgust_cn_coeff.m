function [c_n,c_nC,c_nf,c_nI,alpha_Et,K_f] = BLTgust_cn_coeff(x,gust_states,U,b,beta,alphadot,alphaddot,alpha_tqc,A1,A2,b1,b2,Cg,Cg2,AG,c_n_alpha,f2prime_n,alpha_0L,cn_v,alpha_g,alpha_g0,F_theta0,F_theta0_p,dt)

% Circulatory and inertial gust states
if ~isempty(Cg2)
    gust_states_C = gust_states(1:6);
    gust_states_I = gust_states(7:11);
    CgI = Cg2;
else
    CgI = 0; gust_states_I = 0; gust_states_C = gust_states;
end

% Effective angles of attack: motion-induced, gust-induced and total
alpha_E = beta^2*U/b*(A1*b1*x(1)+A2*b2*x(2))+alpha_tqc*(1-A1-A2);
alpha_gE = Cg*gust_states_C+alpha_g*(1-sum(AG));
alpha_Et = alpha_E+alpha_gE;

% Circulatory, gust-augmented unsteady - attached flow 
c_nC = c_n_alpha*sin(alpha_Et-alpha_0L);

% Circulatory, gust-augmented unsteady - separated flow
K_f = ((1+f2prime_n^(1/2))/2)^2;
c_nf = c_nC*K_f;

% Inertial:  motion-induced, gust-induced and total
c_nIm = pi*b/U^2*(U*alphadot+1/2*b*alphaddot);
c_nIg = CgI*gust_states_I;
if isnan(F_theta0_p), dF_theta0 = 0; else, dF_theta0 = (F_theta0-F_theta0_p); end
c_nIcg = -alpha_g0*b/U*dF_theta0/dt;
c_nI = c_nIm+c_nIg+c_nIcg;

% Total: circulatory + inertial + vortex 
c_n = c_nf+c_nI+cn_v;

end