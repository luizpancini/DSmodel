function [c_m,c_mC,c_mf,c_mI,dCP] = BL_cm_coeff(c_mv,c_nC,c_nf,x,f2prime_m,upstroke,M,U,b,ah,beta,alpha_bar,q,theta,R,RD,S,P,kappa_0,kappa_1,kappa_2,kappa_3,c_m0,K0,K1,K2,A3,b3,A4,b4,b5,K_aM,K_qM,T_I)

% Circulatory unsteady - attached flow (at pitch axis)
c_mC = c_nC*(1/2+ah/2-(1/4-K0))-pi/8*b5*beta^2*U/b*x(7); 

% Center of pressure variables
K1_prime = K1*(1-kappa_1*RD*(1-abs(theta)))-kappa_2*R*abs(theta)*upstroke;  
K2_prime = K2*(1+kappa_3*S*R^2*(1-P)*(~upstroke));                           

% Total center of pressure offset from quarter-chord
dCP = K0+K1_prime*(1-f2prime_m)+K2_prime*sin(pi*f2prime_m^kappa_0);

% Circulatory unsteady - separated flow (at quarter-chord)
c_mf = dCP*c_nf;

% Inertial (at pitch axis)
c_mI = A3*(-2*ah)/(M*b3*K_aM*T_I)*x(5) + A4*(-2*ah)/(M*b4*K_aM*T_I)*x(6) + (ah^2+1/3)/(M*K_qM*T_I)*x(8) - (-2*ah)/M*alpha_bar - (ah^2+1/3)/M*q;

% Total: zero-lift + circulatory + inertial + vortex
c_m = c_m0+c_mf+c_mI+c_mv;

end