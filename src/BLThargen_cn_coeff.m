function [c_n,c_nC,c_nf,c_nI,c_nIa,c_nIf,c_nIg,w_E,w_Ep,w_Ef,w_Eg,alpha_E] = BLThargen_cn_coeff(x,c_nv,gust_states,U,Udot,b,ah,Uc,w_tqc_pp,alpha,alphadot,alphaddot,hddot,deltadot,deltaddot,w_tqc_f,w0,wg,vg,AG,bG,Cg,f2prime_n,T1,T4,alpha_0L,c_n_alpha,eps_fn,F_theta0,F_theta0_i,dt)

% Effective downwash: pitch-plunge-pitch-rate-airspeed-rate-induced, flap-induced, gust-induced, and total
w_Ep = w_tqc_pp-sum(x(8:9))/c_n_alpha;
w_Ef = w_tqc_f-sum(x(10:11));
w_Eg = wg*(1-sum(AG))+U/b*(AG.*bG)'*gust_states;
w_E = (w_Ep+eps_fn*w_Ef+w_Eg);

% Effective circulatory AoA
alpha_E = atan(w_E/(Uc+vg));

% Circulatory unsteady - attached flow
c_nC = c_n_alpha*sin(alpha_E-alpha_0L);

% Circulatory unsteady - separated flow
K_f = ((1+sqrt(f2prime_n))/2)^2;
c_nf = c_nC*K_f;

% Inertial: pitch-plunge-induced, flap-induced, gust-induced and total
c_nIa = pi*b/U^2*(U*alphadot+alpha*Udot+hddot-ah*b*alphaddot);
c_nIf = -eps_fn*b/U^2*(U*T4*deltadot+b*T1*deltaddot);
c_nIg = -w0*b/U^2*(F_theta0-F_theta0_i)/dt;
c_nI = c_nIa+c_nIf+c_nIg;

% Total: circulatory + inertial + vortex
c_n = c_nf+c_nI+c_nv;

end