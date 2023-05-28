function c_n_alpha = get_c_n_alpha_now(airfoil,M)

%% Airfoil parameters table - as functions of Mach number
switch airfoil
    case 'NACA0002'
        % Interpolation strategy
        interp_mode = 'linear';
        % Bound limits
        M = max([0.035,min([0.9,M])]); beta = sqrt(1-M^2);
        % Mach-dependent parameters
        Mach_range =                [0.035;  0.9];
        c_n_alpha_range = 2*pi/beta*[1.00;   1.00]; 
    case 'NACA0006'
        % Interpolation strategy
        interp_mode = 'linear';
        % Bound limits
        M = max([0.035,min([0.9,M])]); beta = sqrt(1-M^2);
        % Mach-dependent parameters
        Mach_range =                [0.035;  0.9];
        c_n_alpha_range = 2*pi/beta*[1.00;   1.00];              
    case 'NACA0012'
        % Interpolation strategy
        interp_mode = 'linear';
        % Bound limits
        M = max([0.035,min([0.3,M])]); beta = sqrt(1-M^2);
        % Mach-dependent parameters
        Mach_range =                [0.035;   0.070;   0.110;   0.185;   0.215;   0.25;    0.28;    0.3];
        c_n_alpha_range = 2*pi/beta*[1.049;   0.972;   0.998;   0.995;   1.00;    1.00;    1.00;    1.00];  
    case 'NACA0018'
        % Interpolation strategy
        interp_mode = 'linear';
        % Bound limits
        M = max([0.060,min([0.160,M])]); beta = sqrt(1-M^2);
        % Mach-dependent parameters
        Mach_range =                [0.062;  0.080;   0.120;   0.150];
        c_n_alpha_range = 2*pi/beta*[0.933;  0.930;   0.941;   0.942];
    case "NACA64A006"
        % Interpolation strategy
        interp_mode = 'linear';
        % Bound limits
        M = max([0.5,min([0.75,M])]); beta = sqrt(1-M^2);
        % Mach-dependent parameters
        Mach_range =                [0.5;   0.75];
        c_n_alpha_range = 2*pi/beta*[1.00;  1.00];       
    otherwise
        error("Airfoil '" + airfoil + "' not listed")
end

%% Interpolated values
% c_n_alpha = interp1(Mach_range,c_n_alpha_range,M,interp_mode);
c_n_alpha = lininterp1(Mach_range,c_n_alpha_range,M);

end