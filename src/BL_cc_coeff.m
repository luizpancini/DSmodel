function c_c = BL_cc_coeff(c_cv,upstroke,alpha,alpha_lag,theta,f2prime_c,alpha_E,R,RD,S,T_s,alpha_0L,alpha1_0c,eta,c_d0,c_n_alpha,E0,E1)

% Lagged AoA to chordwise breakpoint angle ratio
theta_c = abs(alpha_lag)/alpha1_0c;

% Circulatory unsteady  
c_cf = -c_d0*cos(alpha)+(1-eta*RD^2)*c_n_alpha*sin(alpha_E-alpha_0L)^2*f2prime_c^(1/2+theta_c+E0*RD*S*(1-T_s)*R*abs(theta)^(1/2)*~upstroke);
if theta_c > 1
    c_cf = c_cf-E1*(1-RD^3)*min([1,theta_c^3-1]);
end

% Total: circulatory unsteady + vortex-induced 
c_c = c_cf+c_cv;

end