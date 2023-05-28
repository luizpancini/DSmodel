function [K0,K1,K2] = get_CPvars_now(airfoil,M)

%% Airfoil parameters table - as functions of Mach number
switch airfoil
    case 'NACA0002'
        % Interpolation strategy
        interp_mode = 'linear';
        % Mach-dependent parameters
        Mach_range =                [0.01;    0.3;      0.9];
        K0_range =                  [0.00;    0.00;     0.00];   %% Distance between quarter-chord and aerodynamic center (static test) - "hard" parameter
        K1_range =                  [-0.000;  -0.000;   -0.000]; %% Influences movement of center of pressure with separation position (static test) - "hard" parameter
        K2_range =                  [0.00;    0.00;     0.00];   %% Influences movement of center of pressure with separation position (static test) - "hard" parameter
    case 'NACA0006'
        % Interpolation strategy
        interp_mode = 'linear';
        % Bound limits
        M = max([0.035,M]);
        % Mach-dependent parameters
        Mach_range =                [0.035;   0.3;      0.5;      0.65;     0.7;     0.75;    0.8;     0.85;    0.9];
        K0_range =                  [-0.006;  -0.006;   -0.005;   -0.004;   -0.003;  0.000;   0.005;   -0.020;  -0.190];
        K1_range =                  [-0.120;  -0.120;   -0.125;   -0.105;   -0.090;  -0.130;  0.020;   0.020;   0.020]; 
        K2_range =                  [0.04;    0.04;     0.04;     0.04;     0.15;    -0.02;   -0.01;   0.00;    0.00];   
   case 'NACA0012'
        % Interpolation strategy
        interp_mode = 'linear';
        % Bound limits
        M = max([0.035,min([0.3,M])]);
        % Mach-dependent parameters
        Mach_range =                [0.035;   0.070;   0.110;   0.185;   0.215;   0.25;    0.28;    0.3];
        K0_range =                  [0.015;   0.010;   0.000;   0.000;   0.008;   0.010;   0.008;   0.000];
        K1_range =                  [-0.141;  -0.128;  -0.107;  -0.161;  -0.131;  -0.103;  -0.105;  -0.108];
        K2_range =                  [0.010;   0.064;   0.041;   0.064;   0.052;   0.042;   0.040;   0.045]; 
   case {'NACA0012-GU','NACA0015','NACA0015-s','NACA23012A'}
        % Interpolation strategy
        interp_mode = 'linear';
        % Bound limits
        M = max([0.075,min([0.160,M])]);  
       % Mach-dependent parameters
        Mach_range =                [0.078;   0.117;   0.155];
        K0_range =             1e-3*[1.6;     1.0;     1.0];
        K1_range =                  [-0.140;  -0.141;  -0.154];
        K2_range =                  [0.046;   0.049;   0.047];  
   case 'NACA0018'
        % Interpolation strategy
        interp_mode = 'linear';
        % Bound limits
        M = max([0.060,min([0.160,M])]);
        % Mach-dependent parameters
        Mach_range =                [0.062;  0.080;   0.120;   0.150];
        K0_range =                  [0.009;  0.010;   0.008;   0.008];
        K1_range =                  [-0.193; -0.190;  -0.154;  -0.159];
        K2_range =                  [0.063;  0.051;   0.047;   0.047];    
   case "NACA64A006"
        % Interpolation strategy
        interp_mode = 'linear';
        % Bound limits
        M = max([0.5,min([0.75,M])]);
        % Mach-dependent parameters
        Mach_range =                [0.5;       0.75];
        K0_range =                  [-0.010;    -0.010]; 
        K1_range =                  [-0.120;    -0.120]; 
        K2_range =                  [0.02;      0.02];  
    otherwise
        error("Airfoil '" + airfoil + "' not listed")
end

%% Interpolated values
% K0 = interp1(Mach_range,K0_range,M,interp_mode);
% K1 = interp1(Mach_range,K1_range,M,interp_mode);
% K2 = interp1(Mach_range,K2_range,M,interp_mode);
% K3 = interp1(Mach_range,K3_range,M,interp_mode);
K0 = lininterp1(Mach_range,K0_range,M);
K1 = lininterp1(Mach_range,K1_range,M);
K2 = lininterp1(Mach_range,K2_range,M);

end