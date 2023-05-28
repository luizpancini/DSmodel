function [c_m,c_mC,c_mf,c_mI,c_mv,dCP] = BLO_cm_coeff(x,U,b,M,beta,c_nC,c_nf,c_nv,alpha,q,A3,A4,b3,b4,b5,K_aM,K_qM,T_I,K0,K1,K2,kappa,tau_v,TvL,c_m0,x_ac)

% Circulatory  unsteady - attached flow
c_mC = c_nC*(0.25-x_ac)-pi/8*b5*beta^2*U/b*x(7); % Not actually used...

% Impulsive
c_mI = A3/(M*b3*K_aM*T_I)*x(5) + A4/(M*b4*K_aM*T_I)*x(6) + 7/(12*M*K_qM*T_I)*x(8) - 1/M*alpha - 7/(12*M)*q;

% Circulatory unsteady - separated flow
fM = max([x(10); x(12)]);
dCP = (K0 + K1*(1-fM) + K2*sin(pi*fM^kappa));
c_mf = dCP*c_nf; 

% Vortex-induced
CP_v = 0;
if tau_v <= 2*TvL 
    CP_v = 0.25*(1-cos(pi*tau_v/TvL)); 
end
c_mv = -CP_v*c_nv;

% Total
c_m = c_mf+c_mv+c_mI+c_m0;

end
