function [c_n_alpha,K0] = airspeed_vars(M,beta)

% Normal force slope
c_n_alpha = 2*pi/beta;

% Center of pressure movement variables
K0 = CPvars_of_M(M);

end