function [alpha,alphadot,alphaddot,alpha_plunge,alpha_bar,alpha_tqc,q,qdot,r,rdot,qR,qRdot,R,Rdot] = kinematics(t,U,b,ah,k,a_0,a_1,a_1h,k_h,r0)

% Pitch angle and its time derivatives
alpha = a_0+a_1*sin(k*U/b*t);
alphadot = a_1*k*U/b*cos(k*U/b*t);
alphaddot = -a_1*(k*U/b)^2*sin(k*U/b*t);

% Plunge rate 
hdot = U*tan(a_1h)*sin(k_h*U/b*t);

% Composition of static pitch and plunge-induced AoA
alpha_plunge = a_0 + atan(hdot/U);

% Effective pitch-plunge-induced angle of attack at pitch axis and its rates
sa = sin(alpha); ca = cos(alpha);
alpha_bar = atan((U*sa+hdot*ca)/(U*ca-hdot*sa));

% Pitch rate and its time derivative
q = 2*alphadot*b/U;
qdot = 2*alphaddot*b/U;

% Angle of attack at 3/4-chord
alpha_tqc = atan((U*sa+hdot*ca+b*(1/2-ah)*alphadot)/(U*ca-hdot*sa));

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