function params = BLG_airfoil_parameters(params,airfoil,M,U,b)

% Bound limits
if M < 0.035
    M = 0.035;
end

%% Airfoil parameters table - as functions of Mach
if airfoil == "S809"
    Mach_range =                [0.065;  0.10;   0.15];
    alpha_0L_range =     pi/180*[-1.0;   -1.0;   -1.0];
    alpha1_0_range =     pi/180*[8.2;    8.2;    8.2];
    alpha2_0_range =     pi/180*[17.0;   17.0;   17.0];
    delta_alpha1_range = pi/180*[0.0;    0.0;    0.0];
    delta_alpha2_range = pi/180*[0.0;    0.0;    0.0];
    eta_range =                 [1.0;    1.0;    1.0];
    kappa_range =               [6.0;    6.0;    6.0];
    c_d0_range =                [0.012;  0.012;  0.012];
    c_m0_range =                [-0.036; -0.036; -0.036];
    c_n1_range =                [1.90;   1.90;   1.90];
    c_n_alpha_range =    180/pi*[0.112;  0.112;  0.112];
    Df_range =                  [2.0;    2.0;    2.0];
    E0_range =                  [0.00;   0.00;   0.00];
    Ef_range =                  [1.0;    1.0;    1.0];
    f01_range =                 [1.0;    1.0; 	 1.0];
    f02_range =                 [0.0;    0.0; 	 0.0];
    f03_range =                 [0.005;  0.020;	 0.020];
    fb1_range =                 [-0.0217; -0.04995;	 -0.04995];
    fb2_range =                 [3.4934; 2.8844; 2.8844];
    fb3_range =                 [4.97e8; 4.08e4; 4.08e4];
    K0_range =                  [-0.003; -0.003; -0.003];
    K1_range =                  [-0.001; -0.001; -0.001];
    K2_range =                  [-0.025; -0.025; -0.025];
    S1_range =                  [18.269; 12.066; 12.066];
    Tf0_range =                 [3.0;    3.0;    3.0];
    Tp_range =                  [1.7;    1.7;    1.7];
    Tv0_range =                 [6.0;    6.0;    6.0];
    TvL_range =                 [11.0;   11.0;   11.0];
else
    error('airfoil not listed')
end

%% Interpolated values
interp_mode = 'nearest';
params.alpha_0L = interp1(Mach_range,alpha_0L_range,M,interp_mode);
params.alpha1_0 = interp1(Mach_range,alpha1_0_range,M,interp_mode);
params.alpha2_0 = interp1(Mach_range,alpha2_0_range,M,interp_mode);
params.delta_alpha1 = interp1(Mach_range,delta_alpha1_range,M,interp_mode);
params.delta_alpha2 = interp1(Mach_range,delta_alpha2_range,M,interp_mode);
params.eta = interp1(Mach_range,eta_range,M,interp_mode);
params.kappa = interp1(Mach_range,kappa_range,M,interp_mode);
params.c_d0 = interp1(Mach_range,c_d0_range,M,interp_mode);
params.c_m0 = interp1(Mach_range,c_m0_range,M,interp_mode);
params.c_n1 = interp1(Mach_range,c_n1_range,M,interp_mode);
params.c_n_alpha = interp1(Mach_range,c_n_alpha_range,M,interp_mode);
params.Df = interp1(Mach_range,Df_range,M,interp_mode);
params.E0 = interp1(Mach_range,E0_range,M,interp_mode);
params.Ef = interp1(Mach_range,Ef_range,M,interp_mode);
params.f01 = interp1(Mach_range,f01_range,M,'nearest');
params.f02 = interp1(Mach_range,f02_range,M,'nearest');
params.f03 = interp1(Mach_range,f03_range,M,'nearest');
params.fb1 = interp1(Mach_range,fb1_range,M,'nearest');
params.fb2 = interp1(Mach_range,fb2_range,M,'nearest');
params.fb3 = interp1(Mach_range,fb3_range,M,'nearest');
params.K0 = interp1(Mach_range,K0_range,M,interp_mode);
params.K1 = interp1(Mach_range,K1_range,M,interp_mode);
params.K2 = interp1(Mach_range,K2_range,M,interp_mode);
params.S1 = interp1(Mach_range,S1_range,M,'nearest');
params.S2 = log((params.f01+params.fb1*exp(params.S1*params.alpha1_0)-params.f02)/params.fb2)/params.alpha1_0;
params.S3 = log((params.f02+params.fb2*exp(params.S2*params.alpha2_0)-params.f03)/params.fb3)/params.alpha2_0;
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