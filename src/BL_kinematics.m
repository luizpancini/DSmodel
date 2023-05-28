function [b,beta,Uc,w_QS,w_tqc_pp,wdot_tqc_pp,w_tqc_f,wdot_tqc_f,q_QS,r,rdot,qR,R,c_n_alpha,cnawdot_tqc_pp,wg_QS,Fn_theta0,Fm_theta0] = BL_kinematics(initial_condition,t,airfoil,chord,U,Ma,U_ip,U_ip_dot,Un_tqc,Undot_tqc,alphadot,delta_f,deltadot_f,deltaddot_f,Om,Omdot,Th,c_n_alpha,r0,wg,vg,tvec_fd,c_n_alpha_vec_fd,r_vec_fd)

% Semi-chord
b = chord/2;

% Compressibility factor
beta = sqrt(1-Ma^2);

% Chordwise and normal relative wind velocity and acceleration components
Uc = -U(2); 
Un = U(3); 

% Quasi-steady pitch-plunge-induced downwash at attachment point
w_QS = Un;

% Quasi-steady pitch-plunge-induced downwash at 3/4-chord and its rate
w_tqc_pp = Un_tqc;
wdot_tqc_pp = Undot_tqc;

% Rate of the product of c_n_alpha by w_tqc_pp
% c_n_alpha = get_c_n_alpha_now(airfoil,Ma);
% if initial_condition || length(unique([tvec_fd; t])) <= 1
%     c_n_alpha_dot = 0;
% else
%     c_n_alpha_dot = get_backwards_fin_diff(tvec_fd,c_n_alpha_vec_fd,t,c_n_alpha);
% end
c_n_alpha_dot = 0;
cnawdot_tqc_pp = c_n_alpha*wdot_tqc_pp+c_n_alpha_dot*w_tqc_pp;

% Quasi-steady flap-induced downwash at 3/4-chord and its rate
w_tqc_f = Th(10)*U_ip/pi*delta_f+b*Th(11)/(2*pi)*deltadot_f;
wdot_tqc_f = Th(10)/pi*(U_ip*deltadot_f+U_ip_dot*delta_f)+b*Th(11)/(2*pi)*deltaddot_f;

% Effective nondimensional pitch rate
q_QS = alphadot*chord/U_ip;

% Nondimensional effective reduced pitch rate and its rate
r = q_QS/2;
% if initial_condition || length(unique([tvec_fd; t])) <= 1
%     rdot = chord/U_ip^2*(U_ip*Omdot(1)-alphadot*U_ip_dot);
% else
%     rdot = get_backwards_fin_diff(tvec_fd,r_vec_fd,t,r);
% end
rdot = 0;

% Unsigned ratio of reduced pitch rate to critical pitch rate
qR = abs(r)/r0;

% Unsigned capped reduced pitch rate ratio 
R = min([1, qR]);

% Gust kinematics
% if any([vg; wg])
%     [wg_QS,Fn_theta0,Fm_theta0] = gust_kin_inc(t,U_ip,b,wg,vg);
% else
    wg_QS = 0; Fn_theta0 = 0; Fm_theta0 = 0;
% end
end