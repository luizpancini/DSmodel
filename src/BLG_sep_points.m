function [f,fprime,alpha_lag] = BLG_sep_points(c_n_prime,alpha,alpha1_0,alpha1,alpha2_0,alpha2,c_n_alpha,f01,f02,f03,fb1,fb2,fb3,S1,S2,S3)

% Steady condition separation point based on AoA 
if abs(alpha) <= alpha1_0
    f = f01 + fb1*exp(S1*abs(alpha));
elseif abs(alpha) >= alpha2_0
    f = f03 + fb3*exp(S3*abs(alpha));
else
    f = f02 + fb2*exp(S2*abs(alpha));
end

% Lagged AoA
alpha_lag = abs(c_n_prime)/c_n_alpha;

% Unsteady lagged separation point based on lagged AoA
if alpha_lag <= alpha1
    fprime = f01 + fb1*exp(S1*alpha_lag);
elseif alpha_lag >= alpha2
    fprime = f03 + fb3*exp(S3*alpha_lag);
else
    fprime = f02 + fb2*exp(S2*alpha_lag);
end

end