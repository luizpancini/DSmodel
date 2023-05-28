function [U,Udot,M,beta,alpha,alphadot,alphaddot,w_tqc,w_tqc_dot,q,qdot,r,rdot,qR,qRdot,R,Rdot] = tvfs_kinematics(t,U_0,U_1,k_U,b,k,a_0,a_1,r0,a_inf,psi_a)

% Freestream velocity and its rate
U = U_0+U_1*sin(k_U*U_0/b*t);
Udot = U_1*k_U*U_0/b*cos(k_U*U_0/b*t);

% Mach and compressibility factor
M = U/a_inf;
beta = sqrt(1-M^2);

% Angle of attack and its time derivatives
alpha = a_0+a_1*sin(k*U_0/b*t+psi_a); sa = sin(alpha); ca = cos(alpha);
alphadot = a_1*k*U_0/b*cos(k*U_0/b*t+psi_a);
alphaddot = -a_1*(k*U_0/b)^2*sin(k*U_0/b*t+psi_a);

% Pitch rate and its time derivative
q = 2*alphadot*b/U;
qdot = 2*b*(alphaddot/U-alphadot*Udot/U^2);

% Angle of attack at 3/4-chord and its rate
ah = -1/2;
alpha_tqc = atan((U*sa+b*(1/2-ah)*alphadot)/(U*ca));
alphadot_tqc = 1/(tan(alpha_tqc)^2+1)*((Udot*sa+U*alphadot*ca+b*(1/2-ah)*alphaddot)*(U*ca)+(U*sa+b*(1/2-ah)*alphadot)*(U*alphadot*sa-Udot*ca))/(U*ca)^2;

% Downwash at 3/4-chord and its rate
w_tqc = U*alpha_tqc;
w_tqc_dot = U*alphadot_tqc+Udot*alpha_tqc;

% Nondimensional reduced pitch rate and its time derivative
r = q/2;
rdot = qdot/2;

% Unsigned ratio of reduced pitch rate to critical pitch rate and its time derivative
qR = abs(r)/r0; 
qRdot = rdot/(r0*k*U/b);

% Unsigned capped reduced pitch rate ratio and its signed time derivative
R = min([1, qR]);
Rdot = sign(qRdot)*min([1, abs(qRdot)]);

end