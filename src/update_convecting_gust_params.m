function [bG,G,Cg] = update_convecting_gust_params(U,U_0,lambda_g0,b,beta,AG,bG,gust_ind)

% Assuming the gust convection speed is constant, but the freestream is
% not, the gust speed ratio is given as:
lambda_g = U/(U+U_0*(1-lambda_g0)/lambda_g0);

% Set structure with current gust speed ratio
params.lambda_g = lambda_g;

% Get current convecting gust parameters
params = convecting_gust_params(params,gust_ind);
bG = [params.g1; params.g2; bG(3:end)];
G = params.G;
Cg = U/b*beta^2*(AG.*bG(3:end))';

end