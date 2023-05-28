function [params,data] = load_hargen_test(input)

%% Initialize optional test condition variables to default values
% Harmonic pitch oscillation 
a_0 = deg2rad(0);           % Mean angle of pitch angle harmonic oscillation [rad]
a_1 = deg2rad(0);           % Amplitude of pitch angle harmonic oscillation  [rad]
k = 0;                      % Reduced frequency of pitch angle harmonic oscillation []
% Harmonic plunge oscillation 
a_1h = deg2rad(0);          % Amplitude of "plunge angle" harmonic oscillation  [rad]
k_h = 0;                    % Reduced frequency of "plunge angle" harmonic oscillation []
psi_ha = deg2rad(0);        % Phase angle between "plunge angle" and pitch angle harmonic oscillations [rad]
% Harmonically time-varying freestream
update_params = 0;          % Time-step frequency to update airspeed-dependent airfoil parameters (can be any non-negative integer)
lambda_U = 0;               % Normalized amplitude of airspeed harmonic oscillation []
k_U = 0;                    % Reduced frequency of airspeed harmonic oscillation []
psi_Ua = deg2rad(0);        % Phase angle between airspeed and pitch angle harmonic oscillations [rad]
% Harmonic flap oscillation 
dh = 0.5;                   % Semichord-normalized flap hinge position after midchord
d_0 = deg2rad(0);           % Mean angle of flap angle harmonic oscillation [rad]
d_1 = deg2rad(0);           % Amplitude of flap angle harmonic oscillation  [rad]
k_f = 0;                    % Reduced frequency of flap angle harmonic oscillation []
psi_fa = deg2rad(0);        % Phase angle between flap and pitch angle harmonic oscillations [rad]
% Gust
gust_profile = '';          % Gust profile
lambda_g = 1;               % Gust convection speed ratio [] (lambda_g = U/(U+U_g))
plot_increment = true;      % Option to plot coefficient's increments
gust_options.wg_0 = 0;      % Gust normal velocity [m/s]
gust_options.H = 0;         % Normalized gust length [semi-chords]
gust_options.f_sin = 0;     % Gust frequency for sinusoidal gusts [Hz]
gust_options.f_i = 0;       % Initial gust frequency for swept sinusoidal gusts [Hz]
gust_options.f_f = 0;       % Final gust frequency for swept sinusoidal gusts [Hz]
gust_options.tg_end = 0;    % Final gust time for swept sinusoidal gusts [s]
% Misc
time_plot_style = "default";% Style for plots of coefficients vs. time
mod_ind_params = false;     % Option to use different set of circulatory indicial parameters
t_init = 0;                 % Simulation time before gust encounter
tspan = [];                 % Simulation time for gust test
    
%% Set general test conditions for current case
if isnumeric(input)
    switch input
        case 1 % NASA frame 10022
            airfoil = 'NACA0012';
            b = 0.305;
            ah = -1/2;
            a_inf = 340;
            M = 0.301;
            a_0 = deg2rad(12);
            a_1 = deg2rad(9.9);
            k = 0.098;
            time_plot_style = "pitch";
        case 2 % Liiva et al's experiment
            airfoil = 'Vertol-23010';
            b = 6.38*0.0254/2;
            ah = -1/2;
            a_inf = 340;
            M = 0.4;
            a_0 = deg2rad(12.45);
            a_1h = deg2rad(3.14);
            k_h = 0.116;
            time_plot_style = "pitch";
    end
    % Set dependent variables
    beta = sqrt(1-M^2);     % Prandtl-Glauert compressibility factor []
    U = M*a_inf;            % Freestream airspeed [m/s]
    % Set reference data structure
    data.authors = {[]};
else
    % Source of reference data
    source = regexpi(input,'\w*(?=-)','match');
    case_now = str2double(regexpi(input,'(?<=-)\d*','match'));
    switch source
        case "NASA"
            [params,data] = load_NASA(case_now);
            [M,U,b,a_inf,beta,a_0,a_1,k,airfoil,ah] = struct2vars(params);
            time_plot_style = "pitch";
        case "GU"
            [params,data] = load_GU(case_now);
            [M,U,b,a_inf,beta,a_0,a_1,k,airfoil,ah] = struct2vars(params);
            time_plot_style = "pitch";
        case "OSU"
            [params,data] = load_OSU(case_now);
            [M,U,b,a_inf,beta,a_0,a_1,k,airfoil,ah] = struct2vars(params);
            time_plot_style = "pitch";
        case "gust_test"
            [params,data] = load_gust_test(case_now);
            [airfoil,ah,a_0,a_1,a_inf,b,beta,k,M,U,lambda_g,gust_profile,gust_options,t_init,tspan,plot_increment] = struct2vars(params);
            time_plot_style = "gust";
        case "flap_test"
            [params,data] = load_flap_test(case_now);
            [M,U,b,ah,dh,a_inf,beta,a_0,a_1,k,d_0,d_1,k_f,airfoil] = struct2vars(params);
            time_plot_style = "pitch";
        case "tvfs_test"
            [params,data] = load_tvfs_test(case_now);
            [M,U,b,a_inf,beta,a_0,a_1,k,~,~,k_U,psi_a,airfoil,lambda_U,ah] = struct2vars(params);
            psi_Ua = -psi_a;
            time_plot_style = "tvfs";
            mod_ind_params = true;  
            update_params = 0;
        otherwise
            error("Unknown source: " + source);
    end
    % Save data source and current case
    data.source = source; data.case_now = case_now;
end

%% Pack data
% Set additional dependent variables
U_0 = U;            % Initial freestream airspeed [m/s]
U_1 = lambda_U*U;   % Amplitude of airspeed harmonic oscillation []

% Set parameters structure
params = variables2struct(struct(),airfoil,b,ah,a_inf,M,U,beta,a_0,a_1,k,a_1h,k_h,psi_ha,update_params,lambda_U,U_0,U_1,k_U,psi_Ua,dh,d_0,d_1,k_f,psi_fa,gust_profile,lambda_g,t_init,plot_increment,gust_options,tspan,time_plot_style,mod_ind_params);

end