function [X,alpha_gE] = BLTgust_recurrence_states(y,dt,b,U,alpha_g,AG,bG)

% Variables at previous time step
alpha_g_p = y(end-4); U_p = y(end-2); X = y(end-1:end);

% Variables' increments
ds = (U+U_p)/(2*b)*dt;
dalpha_g = alpha_g-alpha_g_p;

% Recurrence solution for "deficiency functions" (states)
X(1) = X(1)*exp(-bG(1)*ds)+AG(1)*dalpha_g*exp(-bG(1)*ds/2);
X(2) = X(2)*exp(-bG(2)*ds)+AG(2)*dalpha_g*exp(-bG(2)*ds/2);

% Effective gust-induced angle of attack
alpha_gE = alpha_g-sum(X);

end