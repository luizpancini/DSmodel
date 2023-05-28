function [alphadot_g,wgdot] = get_alphadot_g(wgdot_fun,t,U,Udot,wg)

% Current gust vertical velocity         
wgdot = wgdot_fun(t);
if isnan(wgdot), wgdot = 0; end

% Gust-induced AoA's rate
alphadot_g = (wgdot/U-wg*Udot/U^2)/(wg^2/U^2+1);

end