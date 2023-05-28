function [c_m,c_mf,c_mI,dCP] = BLTflap_cm_coeff(c_nf,alphadot,alphaddot,delta,deltadot,deltaddot,q,T1,T4,T7,T8,T10,T11,eps_m,K0,K1,K1_f,K2,K3,kappa,c_m0,U,b,dh,f2prime_m,theta,R,S,P,RD,cm_v)

% Circulatory unsteady
if theta*q > 0
    del_CP = -R*abs(theta)*(1-f2prime_m); % Additional term on modeling of the movement of the CP on stalled conditions: proportional to (1-f2prime_m), R*abs(theta) to maintain continuity
    K2_prime = K2;
else
    del_CP = 0; 
    K2_prime = K2*(1+2*R*RD*S*(1-P)); % 2*R because proportional to pitch rate, S to apply only if stall occured, (1-P) to reduce effect at light stall
end
K1 = K1*(1-K1_f*RD^(3/2)*(1-abs(theta))); 
dCP = K0 + K1*(1-f2prime_m) + K2_prime*sin(pi*f2prime_m^kappa) + K3*del_CP;
c_mf = dCP*c_nf;

% Impulsive: airfoil-induced, flap-induced and total
c_mIa = -pi*b/(2*U^2)*(U*alphadot+3/8*b*alphaddot);
c_mIf = -eps_m/(2*U^2)*(U^2*(T4+T10)*delta+U*b*(T1-T8-(dh+1/2)*T4+T11/2)*deltadot+b^2*(T7+T1*(dh+1/2))*deltaddot);
c_mI = c_mIa+c_mIf;

% Total: zero-lift + circulatory + impulsive + vortex
c_m = c_m0+c_mf+c_mI+cm_v;

end
