function [alpha,alphadot,alphaddot,alpha_tqc,delta,deltadot,deltaddot,delta_qs,q,qdot,r,rdot,qR,qRdot,R,Rdot] = flap_kinematics(t,U,b,k,a_0,a_1,r0,k_f,d_0,d_1,T10,T11)

% Angle of attack and its time derivatives
alpha = a_0+a_1*sin(k*U/b*t);
alphadot = a_1*k*U/b*cos(k*U/b*t);
alphaddot = -a_1*(k*U/b)^2*sin(k*U/b*t);

% Flap angle and its time derivatives
delta = d_0+d_1*sin(k_f*U/b*t);
deltadot = d_1*k_f*U/b*cos(k_f*U/b*t);
deltaddot = -d_1*(k_f*U/b)^2*sin(k_f*U/b*t);

% Nondimensional pitch rate and its time derivative
q = 2*alphadot*b/U;
qdot = 2*alphaddot*b/U;

% Angle of attack at 3/4-chord
alpha_tqc = alpha+q/2;

% Quasi-steady effective flap angle
delta_qs = T10/pi*delta+b*T11/(2*pi*U)*deltadot;

% Reduced pitch rate and its time derivative
r = q/2;
rdot = qdot/2;

% Unsigned ratio of reduced pitch rate to critical pitch rate and its time derivative
qR = abs(r)/r0; 
qRdot = rdot/(r0*k*U/b);

% Unsigned capped reduced pitch rate ratio and its signed time derivative
R = min([1, qR]);
Rdot = sign(qRdot)*min([1, abs(qRdot)]);

end