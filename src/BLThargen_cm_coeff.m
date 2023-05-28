function [c_m,c_mC,c_mf,c_mQS,c_mI,c_mIa,c_mIf,c_mIg,dCP] = BLThargen_cm_coeff(airfoil,c_nC,c_nf,c_mv,f2prime_m,theta,R,S,P,RD,upstroke,ah,dh,U,Udot,b,M,alpha,alphadot,alphaddot,hddot,delta,deltadot,deltaddot,Ucm_QS,eps_m,kappa_0,kappa_1,kappa_2,kappa_3,c_m0,T1,T4,T7,T8,T10,T11,w0,F2_theta0,F2_theta0_i,dt)

% Get current value for CP variables
[K0,K1,K2] = get_CPvars_now(airfoil,M);

% Circulatory unsteady - attached flow (at pitch axis)
c_mC = c_nC*(1/2+ah/2-(1/4-K0)); 

% Center of pressure dynamic variables
K1_prime = K1*(1-kappa_1*RD*(1-abs(theta)))-kappa_2*R*abs(theta)*upstroke;  
K2_prime = K2*(1+kappa_3*S*R^2*(1-P)*(~upstroke));                           

% Total center of pressure offset from quarter-chord
dCP = K0+K1_prime*(1-f2prime_m)+K2_prime*sin(pi*f2prime_m^kappa_0);

% Circulatory unsteady - separated flow (at quarter-chord)
c_mf = dCP*c_nf;

% Quasi-steady - at pitch axis
% c_mQS = 0;
c_mQS = Ucm_QS/U;

% Inertial: pitch-plunge-induced, flap-induced, gust-induced and total
c_mIa = -pi*b/(2*U^2)*((1/8+ah^2)*b*alphaddot-ah*hddot); 
c_mIf = -eps_m/(2*U^2)*(U^2*(T4+T10)*delta+U*b*(T1-T8-(dh-ah)*T4+T11/2)*deltadot+b^2*(T7+T1*(dh-ah))*deltaddot); % Includes quasi-steady
c_mIg = w0*b/U^2*(F2_theta0-F2_theta0_i)/dt;
c_mI = c_mIa+c_mIf+c_mIg;

% Total: zero-lift + circulatory + quasi-steady + inertial + vortex 
c_m = c_m0+c_mf+c_mQS+c_mI+c_mv;

end