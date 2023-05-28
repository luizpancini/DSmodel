function [U,Udot,M,Mdot,beta,w_QS,Uc,alpha_QS,alphadot_QS,alpha,alphadot,alphaddot,hdot,hddot,alpha_plunge,delta,deltadot,deltaddot,w_tqc_pp,wdot_tqc_pp,cnawdot_tqc,q,q_QS,r,rdot,qR,R,w_tqc_f,wdot_tqc_f,Ucm_QS,Ucm_QS_dot,c_n_alpha] = ...
         BLhargen_kinematics(t,dt,airfoil,b,ah,a_inf,a_0,a_1,k,a_1h,k_h,psi_ha,U_0,U_1,k_U,psi_Ua,d_0,d_1,k_f,psi_fa,T10,T11,r0,r_i,c_n_alpha_i)
     
% Pitch angle and its rates
alpha = a_0+a_1*sin(k*U_0/b*t); sa = sin(alpha); ca = cos(alpha);
alphadot = a_1*k*U_0/b*cos(k*U_0/b*t);
alphaddot = -a_1*(k*U_0/b)^2*sin(k*U_0/b*t);

% Plunge rate and its rates
hdot = U_0*tan(a_1h)*sin(k_h*U_0/b*t+psi_ha);
hddot = U_0*tan(a_1h)*(k_h*U_0/b)*cos(k_h*U_0/b*t+psi_ha);
% hdddot = -U_0*tan(a_1h)*(k_h*U_0/b)^2*sin(k_h*U_0/b*t+psi_ha);

% Freestream airspeed and its rate
U = U_0+U_1*sin(k_U*U_0/b*t+psi_Ua);
Udot = U_1*k_U*U_0/b*cos(k_U*U_0/b*t+psi_Ua);

% Composition of static pitch and plunge-induced AoA
alpha_plunge = a_0+atan(hdot/U);

% Quasi-steady pitch-plunge-induced downwash and chordwise velocity at pitch axis and their rates
w_QS = U*sa+hdot*ca;
wdot_QS = U*ca*alphadot+Udot*sa+hddot*ca-hdot*sa*alphadot;
Uc = U*ca-hdot*sa;
Ucdot = Udot*ca-U*sa*alphadot-hddot*sa-hdot*ca*alphadot;

% Quasi-steady pitch-plunge-induced angle of attack at pitch axis and its rates
alpha_QS = atan(w_QS/Uc);
alphadot_QS = (wdot_QS*Uc-w_QS*Ucdot)/(w_QS^2+Uc^2);
% alphaddot_QS_approx = alphaddot + (U*hdddot*(U^2+hdot^2)-2*U*hdot*hddot^2)/(U^2+hdot^2)^2; % Approximate, considering steady freestream airspeed and small pitch angle

% Flap angle and its rates
delta = d_0+d_1*sin(k_f*U/b*t+psi_fa);
deltadot = d_1*k_f*U/b*cos(k_f*U/b*t+psi_fa);
deltaddot = -d_1*(k_f*U/b)^2*sin(k_f*U/b*t+psi_fa);

% Mach number, its rate, and compressibility factor
M = U/a_inf;
Mdot = Udot/a_inf;
beta = sqrt(1-M^2);

% Nondimensional pitch rate and its rate
q = 2*alphadot*b/U;
% qdot = 2*b*(alphaddot/U-alphadot*Udot/U^2);

% Effective nondimensional pitch rate and its rate
q_QS = 2*alphadot_QS*b/U;
% qdot_bar = 2*b*(alphaddot_bar_approx/U-alphadot_bar*Udot/U^2);

% Nondimensional effective reduced pitch rate and its rate
r = q_QS/2;
% rdot = qdot_bar/2;
rdot = (r-r_i)/dt;

% Unsigned ratio of reduced pitch rate to critical pitch rate and its reduced rate
qR = abs(r)/r0; 
% qRdot = rdot/(r0*k*U/b);

% Unsigned capped reduced pitch rate ratio and its signed rate
R = min([1, qR]);
% Rdot = sign(qRdot)*min([1, abs(qRdot)]);

% Angle of attack at 3/4-chord and its rate
% alpha_tqc = atan((U*sa+hdot*ca+b*(1/2-ah)*alphadot)/(U*ca-hdot*sa));
% alphadot_tqc = 1/(tan(alpha_tqc)^2+1)*((Udot*sa+U*alphadot*ca+hddot*ca-hdot*alphadot*sa+b*(1/2-ah)*alphaddot)*(U*ca-hdot*sa)+(U*sa+hdot*ca+b*(1/2-ah)*alphadot)*(U*alphadot*sa-Udot*ca+hddot*sa+hdot*alphadot*ca))/(U*ca-hdot*sa)^2;

% Pitch-plunge-pitch-rate-induced downwash at 3/4-chord and its rate
% w_tqc = U*sin(alpha_tqc);
% wdot_tqc = U*cos(alpha_tqc)*alphadot_tqc+Udot*sin(alpha_tqc);
w_tqc_pp = U*sa+hdot*ca+b*(1/2-ah)*alphadot;
wdot_tqc_pp = (U*ca-hdot*sa)*alphadot+hddot*ca+b*(1/2-ah)*alphaddot+Udot*sa;

% Rate of the product of c_n_alpha by w_tqc
if dt>0 && any(Udot), c_n_alpha = get_c_n_alpha_now(airfoil,M); c_n_alpha_dot = (c_n_alpha-c_n_alpha_i)/dt; else, c_n_alpha = c_n_alpha_i; c_n_alpha_dot = 0; end
cnawdot_tqc = c_n_alpha*wdot_tqc_pp+c_n_alpha_dot*w_tqc_pp;

% % Quasi-steady flap-induced angle of attack at 3/4-chord and its rate
% delta_QS = T10/pi*delta+b*T11/(2*pi*U)*deltadot;
% deltadot_QS = T10/pi*deltadot+b*T11/(2*pi)*(deltaddot/U-deltadot*Udot/U^2);

% Quasi-steady flap-induced downwash at 3/4-chord and its rate
w_tqc_f = T10*U/pi*delta+b*T11/(2*pi)*deltadot;
wdot_tqc_f = T10/pi*(U*deltadot+Udot*delta)+b*T11/(2*pi)*deltaddot;

% Product of quasi-steady cm by U and its rate (notice that Fung and Hodges
% & Pierce give the factor as 1/4, whereas Leishman as 1/8. Indeed, from
% the definition of the c_m = M/(2*rho*U^2*b^2) and the center of pressure
% for the corresponding normal force being at the 3/4-chord, the factor
% should be 1/4)
Ucm_QS = -1/4*c_n_alpha*b*(1/2-ah)*alphadot; 
Ucm_QS_dot = -1/4*b*(1/2-ah)*(c_n_alpha*alphaddot+c_n_alpha_dot*alphadot);

end