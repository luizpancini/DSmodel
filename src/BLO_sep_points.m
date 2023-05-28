function [f,fprime,alpha_lag] = BLO_sep_points(so_state,alpha,alpha1_0,alpha1,f0,fb,S1,S2,c_n_alpha)

% Steady condition separation point based on AoA 
if abs(alpha) <= alpha1_0
    f = 1-(1-fb)*exp((abs(alpha)-alpha1_0)/S1);
else
    f = f0+(fb-f0)*exp((alpha1_0-abs(alpha))/S2);
end

% Lagged AoA
alpha_lag = so_state/c_n_alpha;

% Unsteady lagged separation point based on lagged AoA
if alpha_lag <= alpha1
    fprime = 1-(1-fb)*exp((alpha_lag-alpha1)/S1);
else
    fprime = f0+(fb-f0)*exp((alpha1-alpha_lag)/S2);
end

end