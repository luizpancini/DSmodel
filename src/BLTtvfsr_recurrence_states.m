function [X,w_tqc_eff] = BLTtvfsr_recurrence_states(y,t,t_i,b,U,alpha,alphadot,w_tqc,A1,A2,b1,b2)

% Variables at previous time step
alpha_p = y(1); U_p = y(end-3); alphadot_p = y(end-2); X = y(end-1:end);

% Variables' increments
dt = t-t_i;
dalpha = alpha-alpha_p;
dalphadot = alphadot-alphadot_p;
dU = U-U_p;
ds = (U+U_p)/(2*b)*dt;
dw = alpha*dU+U*dalpha+b*dalphadot;

% Recurrence solution for "deficiency functions" (states)
X(1) = X(1)*exp(-b1*ds)+A1*dw*exp(-b1*ds/2);
X(2) = X(2)*exp(-b2*ds)+A2*dw*exp(-b2*ds/2);

% Effective downwash at 3/4-chord
w_tqc_eff = w_tqc-sum(X);

end