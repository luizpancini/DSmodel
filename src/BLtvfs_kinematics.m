function [U,M,Mdot,beta,alpha,alphadot,alphaddot,w_tqc,w_tqc_dot,cnaw_tqc_dot,q,qdot,r,rdot,qR,qRdot,R,Rdot,Ucm_qs,Ucmqs_dot] = ...
         BLtvfs_kinematics(t,U_0,U_1,k_U,b,k,a_0,a_1,r0,a_inf,psi_a)

% Freestream velocity and its rate
U = U_0+U_1*sin(k_U*U_0/b*t);
Udot = U_1*k_U*U_0/b*cos(k_U*U_0/b*t);

% Mach and its rate, compressibility factor and its rate
M = U/a_inf;
Mdot = Udot/a_inf;
beta = sqrt(1-M^2);
betadot = -U/(a_inf^2*beta)*Udot;

% Angle of attack and its rates
alpha = a_0+a_1*sin(k*U_0/b*t+psi_a); sa = sin(alpha); ca = cos(alpha);
alphadot = a_1*k*U_0/b*cos(k*U_0/b*t+psi_a);
alphaddot = -a_1*(k*U_0/b)^2*sin(k*U_0/b*t+psi_a);

% Nondimensional pitch rate and its rate
q = 2*alphadot*b/U;
qdot = 2*b*(alphaddot/U-alphadot*Udot/U^2);

% Angle of attack at 3/4-chord and its rate
alpha_tqc = atan((U*sa+b*alphadot)/(U*ca));
alphadot_tqc = 1/(tan(alpha_tqc)^2+1)*((Udot*sa+U*alphadot*ca+b*alphaddot)*U*ca+(U*sa+b*alphadot)*(U*alphadot*sa-Udot*ca))/(U*ca)^2;

% Downwash at 3/4-chord and its rate
w_tqc = U*alpha_tqc;
w_tqc_dot = U*alphadot_tqc+Udot*alpha_tqc;

% Rate of the product of c_n_alpha = 2*pi/beta by w_tqc
cnaw_tqc_dot = 2*pi/beta*w_tqc_dot-2*pi*betadot/beta^2*w_tqc;

% Nondimensional reduced pitch rate and rate
r = q/2;
rdot = qdot/2;

% Unsigned ratio of reduced pitch rate to critical pitch rate and its rate
qR = abs(r)/r0; 
qRdot = rdot/(r0*k*U/b);

% Unsigned capped reduced pitch rate ratio and its signed time derivative
R = min([1, qR]);
Rdot = sign(qRdot)*min([1, abs(qRdot)]);

% Product of quasi-steady cm by U
Ucm_qs = -1/8*2*pi/beta*b*alphadot;
Ucmqs_dot = -pi/4*b*(alphaddot/beta-alphadot*betadot/beta^2);

end