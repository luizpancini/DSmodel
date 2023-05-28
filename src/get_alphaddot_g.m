function [alphaddot_g,wgddot] = get_alphaddot_g(wgddot_fun,t,U,Udot,wg,wgdot)

% Current gust vertical velocity         
wgddot = wgddot_fun(t);
if isnan(wgddot), wgddot = 0; end

% Gust-induced AoA's rate
alphaddot_g = (U*wgddot+wgdot*Udot)/(U^2+wg^2) - 2*U*wgdot*(wg*wgdot+U*Udot)/(U^2+wg^2)^2;

end