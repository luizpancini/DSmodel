function params = set_cgust_params(params,p)

% Table for parameters as a function of gust convection speed ratio (lambda)
lambda_range = [  0.5;     0.6;    0.7;   0.75;  1.0;     1.5;     2.0];
G_range =      [ p(1);    p(2);   p(3);   p(4);    0;    p(5);    p(6)];
g1_range =     [ p(7);    p(8);   p(9);  p(10);    0;   p(11);   p(12)];
g2_range =     [p(13);   p(14);  p(15);  p(16);    0;   p(17);   p(18)];

% Interpolate
interp_method = 'pchip';
params.G  = interp1(lambda_range,G_range,params.lambda_g,interp_method);
params.g1 = interp1(lambda_range,g1_range,params.lambda_g,interp_method);
params.g2 = interp1(lambda_range,g2_range,params.lambda_g,interp_method);

% Override for better c_n_alpha fit
params.c_n_alpha = 2*pi/params.beta*1.07;

end