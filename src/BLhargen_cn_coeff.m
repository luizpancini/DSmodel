function [c_n,c_nC,c_nf,c_nI,c_nIa,c_nIf,c_nIg,w_E,w_Ep,w_Ef,w_Eg,alpha_E] = BLhargen_cn_coeff(x,gust_states,c_nv,f2prime_n,U,b,beta,Uc,w_tqc_pp,deltadot,deltaddot,w_tqc_f,wg,vg,AG,bG,G,Cg,eps_fn,alpha_0L,c_n_alpha,T1,T4)

% Effective downwash: pitch-plunge-pitch-rate-airspeed-rate-induced, flap-induced, gust-induced, and total
w_Ep = w_tqc_pp-sum(x(8:9))/c_n_alpha; % for alternate states, w_Ep = beta^2*U/b*(A1*b1*x(8)+A2*b2*x(9));
w_Ef = w_tqc_f-sum(x(10:11));
w_Eg = wg*(1-sum(AG))+U/b*beta^2*(AG.*bG(3:end))'*gust_states;
w_E = (w_Ep+eps_fn*w_Ef+w_Eg);

% Effective circulatory AoA
alpha_E = atan(w_E/(Uc+vg));

% Circulatory unsteady - attached flow
c_nC = c_n_alpha*sin(alpha_E-alpha_0L);

% Circulatory unsteady - separated flow
K_f = ((1+sqrt(f2prime_n))/2)^2;
c_nf = c_nC*K_f;

% Inertial: pitch-plunge-airspeed-rate-induced, flap-induced and total
c_nIa = sum(x(12:14));
c_nIf = -eps_fn*b/U^2*(U*T4*deltadot+b*T1*deltaddot);
c_nIg = 1/b*beta^2*G*(bG(1:2).*[1;-1])'*x(21:22);
c_nI = c_nIa+c_nIf+c_nIg;

% Total: circulatory + inertial + vortex
c_n = c_nf+c_nI+c_nv;

end