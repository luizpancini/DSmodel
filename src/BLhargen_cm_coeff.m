function [c_m,c_mC,c_mf,c_mI,c_mQS,dCP] = BLhargen_cm_coeff(airfoil,x,c_nC,c_nf,c_mv,Ucm_QS,f2prime_m,theta,R,S,P,RD,upstroke,ah,dh,U,b,M,delta,deltadot,deltaddot,eps_m,kappa_0,kappa_1,kappa_2,kappa_3,c_m0,T1,T4,T7,T8,T10,T11)

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
% c_mf = c_mC;    % For TVFS tests, c_mf is essentially equal to c_mC
c_mf = dCP*c_nf;

% Quasi-steady - at pitch axis
% c_mQS = 0;
% c_mQS = (Ucm_QS-x(20))/U; % Compressible flow 
c_mQS = Ucm_QS/U;           % Incompressible flow 

% Inertial: pitch-plunge-airspeed-rate-induced, flap-induced and total
c_mIa = sum(x(15:19));
c_mIf = -eps_m/(2*U^2)*(U^2*(T4+T10)*delta+U*b*(T1-T8-(dh-ah)*T4+T11/2)*deltadot+b^2*(T7+T1*(dh-ah))*deltaddot); % Includes quasi-steady
c_mI = c_mIa+c_mIf;

% Total: zero-lift + circulatory + inertial + vortex 
c_m = c_m0+c_mf+c_mQS+c_mI+c_mv;

end