function y = find_f2prime_m(f2prime_m,c_mf_ref,c_nf,R,RD,S,theta,q_bar,kappa_0,kappa_1,kappa_2,kappa_3,K0,K1,K2)

% Center of pressure variables
upstroke = theta*q_bar >= 0;
K1_prime = K1*(1-kappa_1*RD*(1-abs(theta)))-kappa_2*R*abs(theta)*upstroke;  
K2_prime = K2*(1+kappa_3*S*R^2*(~upstroke));                           

% Total center of pressure offset from quarter-chord
dCP = K0 + K1_prime*(1-f2prime_m) + K2_prime*sin(pi*f2prime_m^kappa_0);

% Circulatory unsteady - separated flow (at quarter-chord)
c_mf = dCP*c_nf;

y = c_mf-c_mf_ref;

end