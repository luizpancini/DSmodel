function [c_m,c_mf,c_mI,dCP] = BLtvfsr_cm_coeff(X,c_nf,U,q,Ucm_qs_eff,K0,K1,K1_f,K2,K3,kappa,c_m0,f2prime_m,theta,R,S,P,RD,cm_v)

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

% Quasi-steady
c_mqs = Ucm_qs_eff/U;

% Impulsive
c_mI = sum(X(6:10));

% Total: zero-lift + circulatory + quasi-steady + impulsive + vortex
c_m = c_m0+c_mf+c_mqs+c_mI+cm_v;

end
