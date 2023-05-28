function [alpha_g,alpha,alphadot,alphaddot,alpha_tqc,q,qdot,r,rdot,qR,qRdot,R,Rdot,alpha_bar,q_bar] = gust_kinematics(t,U,b,ah,k,a_0,a_1,r0,wg_fun,wgdot_fun,wgddot_fun)

% Gust-induced AoA and its rates
[alpha_g,wg] = get_alpha_g(wg_fun,t,U);
[alphadot_g,wgdot] = get_alphadot_g(wgdot_fun,t,U,0,wg);
[alphaddot_g,wgddot] = get_alphaddot_g(wgddot_fun,t,U,0,wg,wgdot);

% Pitch angle and its rates
alpha = a_0+a_1*sin(k*U/b*t); 
alphadot = a_1*k*U/b*cos(k*U/b*t);
alphaddot = -a_1*(k*U/b)^2*sin(k*U/b*t);

% Nondimensional pitch rate and its rate
q = 2*alphadot*b/U;
qdot = 2*alphaddot*b/U;

% Effective pitch-gust-induced angle of attack at pitch axis and its rates
ca = cos(alpha); sa = sin(alpha);
alpha_bar = atan((U*sa+wg*ca)/(U*ca-wg*sa));
alphadot_bar = 1/(tan(alpha_bar)^2+1)*((U*alphadot*ca+wgdot*ca-wg*alphadot*sa)*(U*ca-wg*sa)+(U*sa+wg*ca)*(U*alphadot*sa+wgdot*sa+wg*alphadot*ca))/(U*ca-wg*sa)^2;
alphaddot_bar_approx = alphaddot + (U*wgddot*(U^2+wg^2)-2*U*wg*wgdot^2)/(U^2+wg^2)^2; % Approximate, considering small pitch angle (sum individual contributions from pitch and gust)

% Effective nondimensional pitch rate and its rate
q_bar = 2*alphadot_bar*b/U;
qdot_bar = 2*b*alphaddot_bar_approx/U;

%%%%%%%%% If we treat the effective AoA and pitch rate as the
%%%%%%%%% pitch-induced, the correlation for the increment in lift is
%%%%%%%%% better.
% alpha_bar = alpha; q_bar = q; qdot_bar = qdot;

% Pitch-induced angle of attack at 3/4-chord
alpha_tqc = atan((U*sa+b*(1/2-ah)*alphadot)/(U*ca));

% Effective nondimensional reduced pitch rate and its rate
r = q_bar/2;
rdot = qdot_bar/2;

% Unsigned ratio of reduced pitch rate to critical pitch rate and its reduced rate
qR = abs(r)/r0; 
qRdot = rdot/(r0*k*U/b);

% Unsigned capped reduced pitch rate ratio and its signed time derivative
R = min([1, qR]);
Rdot = sign(qRdot)*min([1, abs(qRdot)]);

end