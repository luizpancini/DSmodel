function [c_m,c_mf,c_mI,dCP] = BLTgust_cm_coeff(c_nf,alphadot,alphaddot,q,K0,K1,K1_f,K2,K3,kappa,c_m0,U,b,f2prime_m,theta,R,S,P,RD,cm_v,alpha_g0,F2_theta0,F2_theta0_p,dt)

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

% Inertial: motion-induced, gust-induced and total 
c_mIm = -pi*b/(2*U^2)*(U*alphadot+3/8*b*alphaddot);
if isnan(F2_theta0_p), dF2_theta0 = 0; else, dF2_theta0 = (F2_theta0-F2_theta0_p); end
c_mIg = alpha_g0*b/U*dF2_theta0/dt;
c_mI = c_mIm+c_mIg;

% Total: zero-lift + circulatory + inertial + vortex
c_m = c_m0+c_mf+c_mI+cm_v;

end