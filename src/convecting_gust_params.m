function params = convecting_gust_params(params,gust_ind)

% Table for parameters as a function of gust convection speed ratio (lambda)
switch gust_ind
    case "K"
        lambda_range =  [  0.5;      0.6;       0.7;     0.75;      1.0;      1.5;      2.0];
        G_range =       [ 8.49;     7.28;      9.64;     4.80;        0;    -2.48;    -4.78];
        g1_range =      [3.361;    3.571;     2.871;    2.956;        0;    1.426;    1.154];
        g2_range =      [1.762;    1.513;     2.039;    1.830;        0;    0.329;    0.363];
    case "CFD"
        lambda_range =  [  0.5;      0.6;       0.7;     0.75;      1.0;      1.5;      2.0];
        G_range =       [ 6.70;     6.00;      3.36;     2.82;        0;    -1.90;    -3.90];
        g1_range =      [2.269;    3.359;     3.110;    2.714;        0;    1.900;    1.346];
        g2_range =      [1.477;    1.224;     1.262;    1.278;        0;    0.337;    0.352];
    case {"BR-C","BR-F"}
        lambda_range =  [  0.5;      0.6;       0.7;     0.75;      1.0;      1.5;      2.0];
        G_range =       [10.56;    8.700;     1.017;    1.663;        0;    -2.40;    -5.10];
        g1_range =      [2.712;    2.623;     4.933;    2.573;        0;    2.570;    1.385];
        g2_range =      [1.943;    1.610;     1.303;    1.829;        0;    0.330;    0.387];
end

% Unpack gust convection speed ratio and limit to within interpolation bounds
lambda_g = max([lambda_range(1),min([lambda_range(end),params.lambda_g])]);
if lambda_g < params.lambda_g || lambda_g > params.lambda_g
    warning('Gust convection speed ratio for current case was limited to within bounds for interpolation of coefficients')
end

% Interpolate
interp_method = 'pchip';
params.G  = interp1(lambda_range,G_range,lambda_g,interp_method);
params.g1 = interp1(lambda_range,g1_range,lambda_g,interp_method);
params.g2 = interp1(lambda_range,g2_range,lambda_g,interp_method);

end