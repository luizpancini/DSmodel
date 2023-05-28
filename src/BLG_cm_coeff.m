function [c_m,c_mC,c_mf,c_mI,c_mv,xCP] = BLG_cm_coeff(x,U,b,M,beta,c_nC,c_nf,c_nv,alpha,q,A3,A4,b3,b4,b5,K_aM,K_qM,T_I,K0,K1,K2,kappa,tau_v,TvL,c_m0,x_ac,alpha_lag,alpha1_0,alpha2_0,alpha2,f,fprime,f2prime,f02,fb2,S2)

% Circulatory
c_mC = c_nC*(0.25-x_ac)-pi/8*b5*beta^2*U/b*x(7); % Not actually used...

% Impulsive
c_mI = A3/(M*b3*K_aM*T_I)*x(5) + A4/(M*b4*K_aM*T_I)*x(6) + 7/(12*M*K_qM*T_I)*x(8) - 1/M*alpha - 7/(12*M)*q;

% CP displacement
fM = f2prime;
xCP = K0 + K1*(1-fM) + K2*sin(pi*fM^kappa);
if alpha_lag > alpha2
%     xCP = K0 + K1*exp(K2*fM^kappa); % This is certainly wrong
      xCP = xCP + (alpha_lag-alpha2)*5*K1*exp(3*fM^(-1/kappa)); % (alpha_lag-alpha2) to maintain continuity
end
xCP = max([xCP -0.5]); % Limit CP position to 0.5 

% Circulatory unsteady
c_mf = xCP*c_nf; 

% Vortex-induced
CP_v = 0;
if tau_v <= 2*TvL 
    CP_v = 0.20*(1-cos(pi*tau_v/TvL)); 
end
c_mv = -CP_v*c_nv;

% Total
c_m = c_mf+c_mv+c_mI+c_m0;

end
