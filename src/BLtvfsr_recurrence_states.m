function [X,cnaw_tqc,cnaw_tqc_eff,Ucm_qs_eff] = BLtvfsr_recurrence_states(y,t,t_i,a_inf,b,U,M,beta,alpha,w_tqc,q,Ucm_qs,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,b5)

% Current inertial indicial parameters
[T_a,T_q,T_M,T_am,T_qm,T_Mm] = get_inertial_ind_params(A1,A2,A3,A4,b1,b2,b3,b4,b5,M,beta);

% Current circulatory variables
cnaw_tqc = c_n_alpha*w_tqc;

% Variables at previous time step
alpha_p = y(1); q_p = y(25); U_p = y(36); cnaw_tqc_p = y(38); Ucm_qs_p = y(39); X = y(40:end);

% Variables' increments
dt = t-t_i;
dalpha = alpha-alpha_p;
dU = U-U_p;
dq = q-q_p;
dM = dU/a_inf;
ds = (U+U_p)/(2*b)*dt;
dcnaw_tqc = cnaw_tqc-cnaw_tqc_p;
dUcmqs = Ucm_qs-Ucm_qs_p;

% Recurrence solution for "deficiency functions" (states)
X(1) = X(1)*exp(-b1*beta^2*ds)+A1*dcnaw_tqc*exp(-b1*beta^2*ds/2);
X(2) = X(2)*exp(-b2*beta^2*ds)+A2*dcnaw_tqc*exp(-b2*beta^2*ds/2);
X(3) = X(3)*exp(-ds/T_a)+4/M*dalpha*exp(-ds/2/T_a);
X(4) = X(4)*exp(-ds/T_q)+1/M*dq*exp(-ds/2/T_q);
X(5) = X(5)*exp(-ds/T_M)+4*alpha/M^2*dM*exp(-ds/2/T_M);
X(6) = X(6)*exp(-ds/(b3*T_am))-A3*1/M*dalpha*exp(-ds/2/(b3*T_am));
X(7) = X(7)*exp(-ds/(b4*T_am))-A4*1/M*dalpha*exp(-ds/2/(b4*T_am));
X(8) = X(8)*exp(-ds/(b3*T_Mm))-A3*alpha/M^2*dM*exp(-ds/2/(b3*T_Mm));
X(9) = X(9)*exp(-ds/(b4*T_Mm))-A4*alpha/M^2*dM*exp(-ds/2/(b4*T_Mm));
X(10) = X(10)*exp(-ds/T_qm)-7/(12*M)*dq*exp(-ds/2/T_qm);
X(11) = X(11)*exp(-b5*beta^2*ds)+dUcmqs*exp(-b5*beta^2*ds/2);

% Effective circulatory variables
cnaw_tqc_eff = cnaw_tqc-sum(X(1:2));
Ucm_qs_eff = Ucm_qs-X(11);

end