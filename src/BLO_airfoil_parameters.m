function params = BLO_airfoil_parameters(params,airfoil,M,U,b)

% Bound limits
if M < 0.035
    M = 0.035;
end

%% Airfoil parameters table - as functions of Mach
if strcmp(airfoil,'NACA0012')
    Mach_range =                [0.035;   0.072;   0.110;   0.185;   0.215;   0.25;    0.3;     0.4;     0.5;     0.6;     0.7;     0.75;    0.8];
    alpha_0L_range =     pi/180*[-0.0;    -0.0;    -0.0;    -0.0;    -0.0;    -0.0;    -0.0;    -0.0;    -0.0;    -0.0;    -0.0;    -0.0;    -0.0];
    alpha1_0_range =     pi/180*[12.4;    13.8;    15.2;    16.3;    16.5;    16.0;    13.7;    12.5;    10.5;    8.5;     5.6;     3.5;     0.7];
    delta_alpha1_range = pi/180*[2.0;     1.0;     2.5;     2.5;     2.5;     2.5;     0.5;     2.0;     1.45;    1.0;     0.8;     0.2;     0.1];
    c_d0_range =                [0.010;   0.010;   0.010;   0.010;   0.010;   0.010;   0.010;   0.008;   0.0077;  0.0078;  0.0078;  0.0079;  0.0114];
    c_m0_range =                [-0.0037; -0.0037; -0.0037; -0.0037; -0.0037; -0.0037; -0.0037; -0.0037; -0.0037; -0.0037; -0.0037; -0.0037; -0.0037];
    eta_range =                 [0.95;    0.95;    0.95;    0.95;    0.95;    0.95;    0.95;    0.95;    0.95;    0.95;    0.95;    0.95;    0.95];
    kappa_range =               [2.0;     2.0;     2.0;     2.0;     2.0;     2.0;     2.0;     2.0;     2.0;     2.0;     2.0;     2.0;     2.0];
    c_n1_range =                [1.25;    1.45;    1.60;    1.80;    1.80;    1.85;    1.60;    1.2;     1.05;    0.92;    0.68;    0.5;     0.18];
    c_n_alpha_range =    180/pi*[0.105;   0.108;   0.108;   0.110;   0.113;   0.115;   0.116;   0.113;   0.117;   0.127;   0.154;   0.175;   0.216];
    Df_range =                  [8.0;     8.0;     8.0;     8.0;     8.0;     8.0;     8.0;     7.75;    6.2;     6.0;     5.9;     5.5;     4.0];
    E0_range =                  [0.00;    0.00;    0.00;    0.00;    0.00;    0.00;    0.00;    0.00;    0.00;    0.00;    0.00;    0.00;    0.00];
    f0_range =                  [0.02;    0.02;    0.02;    0.02;    0.02;    0.02;    0.02;    0.02;    0.02;    0.02;    0.02;    0.02;    0.02];
    fb_range =                  [0.70;    0.70;    0.70;    0.70;    0.70;    0.70;    0.70;    0.70;    0.70;    0.70;    0.70;    0.70;    0.70];
    K0_range =                  [0.0025;  0.0025;  0.0025;  0.0025;  0.0025;  0.0025;  0.0025;  0.006;   0.02;    0.038;   0.03;    0.001;   -0.01];
    K1_range =                  [-0.120;  -0.120;  -0.120;  -0.120;  -0.120;  -0.120;  -0.120;  -0.135;  -0.125;  -0.12;   -0.09;   -0.13;   0.02];
    K2_range =                  [0.04;    0.04;    0.04;    0.04;    0.04;    0.04;    0.04;    0.05;    0.04;    0.04;    0.15;    -0.02;   -0.01];
    S1_range =           pi/180*[3.0;     3.0;     3.0;     3.5;     3.5;     3.5;     3.0;     3.25;    3.5;     4.0;     4.5;     3.5;     0.70];
    S2_range =           pi/180*[1.5;     1.5;     1.5;     2.0;     2.0;     2.0;     1.5;     1.6;     1.2;     0.7;     0.5;     0.8;     0.18];
    Tf0_range =                 [3.0;     3.0;     3.0;     3.0;     3.0;     3.0;     3.0;     2.5;     2.2;     2.0;     2.0;     2.0;     2.0];
    Tp_range =                  [1.7;     1.7;     1.7;     1.7;     1.7;     1.7;     1.7;     1.8;     2.0;     2.5;     3.0;     3.3;     4.3];
    Tv0_range =                 [6.0;     6.0;     6.0;     6.0;     6.0;     6.0;     6.0;     6.0;     6.0;     6.0;     6.0;     6.0;     4.0];
    TvL_range =                 [4.0;     5.0;     5.0;     5.0;     5.0;     5.0;     5.0;     9.0;     9.0;     9.0;     9.0;     9.0;     9.0];
elseif strcmp(airfoil,'AMES-01')
    Mach_range =                [0.035;   0.072;   0.110;   0.185;   0.215;   0.25;    0.3];
    alpha_0L_range =     pi/180*[-0.8;    -0.8;    -0.8;    -0.8;    -0.8;    -0.8;    -0.8];
    alpha1_0_range =     pi/180*[12.9;    14.5;    16.0;    16.0;    16.0;    15.3;    14.9];
    delta_alpha1_range = pi/180*[1.0;     1.0;     1.0;     1.0;     1.0;     1.0;     1.0];
    eta_range =                 [0.95;    0.95;    0.95;    0.95;    0.95;    0.95;    0.95];
    kappa_range =               [2.0;     2.0;     2.0;     2.0;     2.0;     2.0;     2.0];
    c_d0_range =                [0.005;   0.005;   0.005;   0.005;   0.005;   0.005;   0.005];
    c_m0_range =                [-0.000;  -0.015;  -0.005;  -0.005;  -0.005;  -0.005;  -0.005];
    c_n1_range =                [1.80;    1.80;    1.80;    1.80;    1.80;    1.80;    1.80];
    c_n_alpha_range =    180/pi*[0.110;   0.110;   0.110;   0.115;   0.118;   0.120;   0.120];
    Df_range =                  [4.0;     4.0;     4.0;     4.0;     4.0;     4.0;     4.0];
    E0_range =                  [0.10;    0.10;    0.10;    0.05;    0.00;   -0.05;   -0.10];
    f0_range =                  [0.01;    0.01;    0.01;    0.01;    0.01;    0.01;    0.01];
    fb_range =                  [0.75;    0.75;    0.75;    0.75;    0.75;    0.75;    0.75];
    K0_range =                  [0.02;    0.02;    0.00;    0.00;    0.00;    0.00;    0.00];
    K1_range =                  [-0.100;  -0.080;  -0.08;   -0.08;   -0.08;   -0.08;   -0.08];
    K2_range =                  [0.02;    0.02;    0.04;    0.04;    0.04;    0.03;    0.03];
    S1_range =           pi/180*[3.0;     3.0;     1.5;     1.5;     1.5;     1.2;     1.2];
    S2_range =           pi/180*[1.5;     1.5;     2.0;     2.0;     2.0;     2.0;     2.0];
    Tf0_range =                 [2.0;     3.0;     2.5;     2.5;     2.5;     2.2;     2.2];
    Tp_range =                  [2.0;     3.0;     2.5;     2.5;     2.5;     2.2;     2.2];
    Tv0_range =                 [2.0;     3.0;     2.5;     2.5;     2.5;     2.2;     2.2];
    TvL_range =                 [4.0;     5.0;     4.5;     4.5;     4.5;     4.2;     4.2];
elseif airfoil == "S809"
    Mach_range =                [0.07;   0.10;   0.15];
    alpha_0L_range =     pi/180*[-0.51;  -0.51;  -1.51];
    alpha1_0_range =     pi/180*[18.0;   19.0;   20.0];
    delta_alpha1_range = pi/180*[1.0;    1.0;    1.0];
    eta_range =                 [1.0;    1.0;    1.0];
    kappa_range =               [2.0;    2.0;    2.0];
    c_d0_range =                [0.012;  0.012;  0.012];
    c_m0_range =                [-0.036; -0.036; -0.036];
    c_n1_range =                [1.90;   1.90;   1.90];
    c_n_alpha_range =    180/pi*[0.112;  0.112;  0.112];
    Df_range =                  [2.0;    2.0;    2.0];
    E0_range =                  [0.10;   0.10;   0.10];
    f0_range =                  [0.01;   0.01; 	 0.01];
    fb_range =                  [0.30;   0.30;	 0.30];
    K0_range =                  [-0.003; -0.003; -0.003];
    K1_range =                  [-0.100; -0.100; -0.100];
    K2_range =                  [-0.025; -0.025; -0.025];
    S1_range =           pi/180*[7.5;    7.5;    7.5];
    S2_range =           pi/180*[2.0;    2.0;    2.0];
    Tf0_range =                 [3.0;    3.0;    3.0];
    Tp_range =                  [1.7;    1.7;    1.7];
    Tv0_range =                 [6.0;    6.0;    6.0];
    TvL_range =                 [11.0;   11.0;   11.0];
else
    error('airfoil not listed')
end

%% Interpolated values
interp_mode = 'linear';
params.alpha_0L = interp1(Mach_range,alpha_0L_range,M,interp_mode);
params.alpha1_0 = interp1(Mach_range,alpha1_0_range,M,interp_mode);
params.delta_alpha1 = interp1(Mach_range,delta_alpha1_range,M,interp_mode);
params.eta = interp1(Mach_range,eta_range,M,interp_mode);
params.kappa = interp1(Mach_range,kappa_range,M,interp_mode);
params.c_d0 = interp1(Mach_range,c_d0_range,M,interp_mode);
params.c_m0 = interp1(Mach_range,c_m0_range,M,interp_mode);
params.c_n1 = interp1(Mach_range,c_n1_range,M,interp_mode);
params.c_n_alpha = interp1(Mach_range,c_n_alpha_range,M,interp_mode);
params.Df = interp1(Mach_range,Df_range,M,interp_mode);
params.E0 = interp1(Mach_range,E0_range,M,interp_mode);
params.f0 = interp1(Mach_range,f0_range,M,interp_mode);
params.fb = interp1(Mach_range,fb_range,M,interp_mode);
params.K0 = interp1(Mach_range,K0_range,M,interp_mode);
params.K1 = interp1(Mach_range,K1_range,M,interp_mode);
params.K2 = interp1(Mach_range,K2_range,M,interp_mode);
params.S1 = interp1(Mach_range,S1_range,M,interp_mode);
params.S2 = interp1(Mach_range,S2_range,M,interp_mode);
params.Tf0 = interp1(Mach_range,Tf0_range,M,interp_mode);
params.Tv0 = interp1(Mach_range,Tv0_range,M,interp_mode);
params.Tp = interp1(Mach_range,Tp_range,M,interp_mode);
params.TvL = interp1(Mach_range,TvL_range,M,interp_mode);
params.x_ac = 0.25-params.K0;

%% Time delay constants adjustment for dimensional time
params.Tf0 = params.Tf0*b/U;
params.Tp = params.Tp*b/U;
params.Tv0 = params.Tv0*b/U;
params.TvL = params.TvL*b/U;

end