function [alpha_g,wg] = get_alpha_g(wg_fun,t,U)

% Current gust vertical velocity       
wg = wg_fun(t);

% Gust-induced angle of attack
alpha_g = atan(wg/U);

end