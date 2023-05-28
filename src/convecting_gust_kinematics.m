function [alpha_g,alpha_g0,alpha,alphadot,alphaddot,alpha_tqc,q,qdot,r,rdot,qR,qRdot,R,Rdot,F_theta0,F2_theta0,theta0] = convecting_gust_kinematics(t,U,b,k,a_0,a_1,r0,wg_fun,wgdot_fun,lambda_g)

% Gust-induced AoA 
[alpha_g0,alpha_g,F_theta0,F2_theta0,theta0] = get_convecting_alpha_g(wg_fun,t,U,b,lambda_g);

% Motion-induced AoA and its rates
alpha = a_0 + a_1*sin(k*U/b*t);
alphadot = a_1*k*U/b*cos(k*U/b*t);
alphaddot = -a_1*(k*U/b)^2*sin(k*U/b*t);

% Motion-induced nondimensional pitch rate and its time derivative
q = 2*alphadot*b/U;
qdot = 2*alphaddot*b/U;

% Angle of attack at 3/4-chord
alpha_tqc = alpha+q/2;

% Motion-induced nondimensional reduced pitch rate and its time derivative
r = q/2;
rdot = qdot/2;

% Unsigned ratio of reduced pitch rate to critical pitch rate and its time derivative
qR = abs(r)/r0; 
qRdot = rdot/(r0*k*U/b);

% Unsigned capped reduced pitch rate ratio and its signed time derivative
R = min([1, qR]);
Rdot = sign(qRdot)*min([1, abs(qRdot)]);

end