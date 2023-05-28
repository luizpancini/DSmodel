clc
clear
% close all

%% Select case
case_now = 6;

%% Load experimental conditions
switch case_now
    case 1
        % Experimental data
        filename_cn = 'CN_Vertol23010_M_04_k_0068_a0_125_a1h_25_exp_raw.mat';
        filename_cm = 'CM_Vertol23010_M_04_k_0068_a0_125_a1h_25_exp_raw.mat';
        load(filename_cn); load(filename_cm);
        dh_cn_exp = CN_Vertol23010_M_04_k_0068_a0_125_a1h_25_exp_new(:,1);
        dh_cm_exp = CM_Vertol23010_M_04_k_0068_a0_125_a1h_25_exp_new(:,1);
        cn_exp = CN_Vertol23010_M_04_k_0068_a0_125_a1h_25_exp_new(:,2);
        cm_exp = CM_Vertol23010_M_04_k_0068_a0_125_a1h_25_exp_new(:,2);
        % Test conditions
        M = 0.4;     % Mach number
        f = 17.24;   % Dimensional frequency [Hz]
        a0 = 12.25;  % Pitch angle [deg]
    case 2
        % Experimental data
        filename_cn = 'CN_Vertol23010_M_04_k_0068_a0_145_a1h_25_exp_raw.mat';
        filename_cm = 'CM_Vertol23010_M_04_k_0068_a0_145_a1h_25_exp_raw.mat';
        load(filename_cn); load(filename_cm);
        dh_cn_exp = CN_Vertol23010_M_04_k_0068_a0_145_a1h_25_exp_raw(:,1);
        dh_cm_exp = CM_Vertol23010_M_04_k_0068_a0_145_a1h_25_exp_raw(:,1);
        cn_exp = CN_Vertol23010_M_04_k_0068_a0_145_a1h_25_exp_raw(:,2);
        cm_exp = CM_Vertol23010_M_04_k_0068_a0_145_a1h_25_exp_raw(:,2);
        % Test conditions
        M = 0.4;     % Mach number
        f = 17.24;   % Dimensional frequency [Hz]
        a0 = 14.65;  % Pitch angle [deg]
    case 3
        % Experimental data
        filename_cn = 'CN_Vertol23010_M_04_k_0121_a0_125_a1h_25_exp_raw.mat';
        filename_cm = 'CM_Vertol23010_M_04_k_0121_a0_125_a1h_25_exp_raw.mat';
        load(filename_cn); load(filename_cm);
        dh_cn_exp = CN_Vertol23010_M_04_k_0121_a0_125_a1h_25_exp_raw(:,1);
        dh_cm_exp = CM_Vertol23010_M_04_k_0121_a0_125_a1h_25_exp_raw(:,1);
        cn_exp = CN_Vertol23010_M_04_k_0121_a0_125_a1h_25_exp_raw(:,2);
        cm_exp = CM_Vertol23010_M_04_k_0121_a0_125_a1h_25_exp_raw(:,2);
        % Test conditions
        M = 0.4;     % Mach number
        f = 30.96;   % Dimensional frequency [Hz]
        a0 = 12.46;  % Pitch angle [deg]
        % From Tyler and Leishman
        load('CN_Vertol23010_M_04_k_0116_a0_12_a1h_3_exp.mat');
        alpha_Tyler = CN_Vertol23010_M_04_k_0116_a0_12_a1h_3_exp(:,1);
        cn_Tyler = CN_Vertol23010_M_04_k_0116_a0_12_a1h_3_exp(:,2);
    case 4
        % Experimental data
        filename_cn = 'CN_Vertol23010_M_04_k_0128_a0_15_a1h_25_exp_raw.mat';
        filename_cm = 'CM_Vertol23010_M_04_k_0128_a0_15_a1h_25_exp_raw.mat';
        load(filename_cn); load(filename_cm);
        dh_cn_exp = CN_Vertol23010_M_04_k_0128_a0_15_a1h_25_exp_raw(:,1);
        dh_cm_exp = CM_Vertol23010_M_04_k_0128_a0_15_a1h_25_exp_raw(:,1);
        cn_exp = CN_Vertol23010_M_04_k_0128_a0_15_a1h_25_exp_raw(:,2);
        cm_exp = CM_Vertol23010_M_04_k_0128_a0_15_a1h_25_exp_raw(:,2);
        % Test conditions
        M = 0.4;     % Mach number
        f = 32.89;   % Dimensional frequency [Hz]
        a0 = 14.88;  % Pitch angle [deg]
    case 5
        % Experimental data
        filename_cn = 'CN_NACA0012_M_02_k_024_a0_125_plunge_exp_raw.mat';
        filename_cm = 'CM_NACA0012_M_02_k_024_a0_125_plunge_exp_raw.mat';
        load(filename_cn); load(filename_cm);
        dh_cn_exp = CN_NACA0012_M_02_k_024_a0_125_plunge_exp_raw(:,1);
        dh_cm_exp = CM_NACA0012_M_02_k_024_a0_125_plunge_exp_raw(:,1);
        cn_exp = CN_NACA0012_M_02_k_024_a0_125_plunge_exp_raw(:,2);
        cm_exp = CM_NACA0012_M_02_k_024_a0_125_plunge_exp_raw(:,2);
%         cn_exp = flipud(cn_exp);
        cm_exp = flipud(cm_exp);
        % Test conditions
        M = 0.2;     % Mach number
        f = 32.47;   % Dimensional frequency [Hz]
        a0 = 12.36;  % Pitch angle [deg]
    case 6
        % Experimental data
        filename_cn = 'CN_NACA0012_M_02_k_024_a0_15_plunge_exp_raw.mat';
        filename_cm = 'CM_NACA0012_M_02_k_024_a0_15_plunge_exp_raw.mat';
        load(filename_cn); load(filename_cm);
        dh_cn_exp = CN_NACA0012_M_02_k_024_a0_15_plunge_exp_raw(:,1);
        dh_cm_exp = CM_NACA0012_M_02_k_024_a0_15_plunge_exp_raw(:,1);
        cn_exp = CN_NACA0012_M_02_k_024_a0_15_plunge_exp_raw(:,2);
        cm_exp = CM_NACA0012_M_02_k_024_a0_15_plunge_exp_raw(:,2);
%         cn_exp = flipud(cn_exp);
%         cm_exp = flipud(cm_exp);
        % Test conditions
        M = 0.2;     % Mach number
        f = 31.75;   % Dimensional frequency [Hz]
        a0 = 14.67;  % Pitch angle [deg]    
    case 8
        % load('CN_Vertol23010_M_04_k_0058_a0_13_a1h_2_exp.mat');
        % alpha_Tyler = CN_Vertol23010_M_04_k_0058_a0_13_a1h_2_exp(:,1);
        % cn_Tyler = CN_Vertol23010_M_04_k_0058_a0_13_a1h_2_exp(:,2);
        % load('plunge_cn3.mat');
        % dh = plunge_cn3(:,1);
        % cn = plunge_cn3(:,2);
        % % Experimental conditions
        % a0 = 12.53; f = 30.96; M = 0.4;
    case 9
        % load('CN_Vertol23010_M_04_k_0058_a0_13_a1h_2_exp.mat');
        % alpha_Tyler = CN_Vertol23010_M_04_k_0058_a0_13_a1h_2_exp(:,1);
        % cn_Tyler = CN_Vertol23010_M_04_k_0058_a0_13_a1h_2_exp(:,2);
        % load('plunge_cn2.mat');
        % dh = plunge_cn2(:,1);
        % cn = plunge_cn2(:,2);
        % % Experimental conditions
        % a0 = 12.53; f = 17.24; M = 0.4;
    case 10
        % load('CN_Vertol23010_M_04_k_0116_a0_12_a1h_3_exp.mat');
        % alpha_Tyler = CN_Vertol23010_M_04_k_0116_a0_12_a1h_3_exp(:,1);
        % cn_Tyler = CN_Vertol23010_M_04_k_0116_a0_12_a1h_3_exp(:,2);
        % load('plunge_cn2.mat');
        % dh = plunge_cn2(:,1);
        % cn = plunge_cn2(:,2);
end

%% Additional test conditions
a_inf = 325;                                                        % Sound speed [m/s]
b = 6.38*0.0254/2;                                                  % Semichord [m]
U = M*a_inf;                                                        % Airspeed [m/s]
T = 1/f;                                                            % Period [s]
k_h = 2*pi*f*b/U;                                                   % Reduced frequency [] 
a1_h = atand(k_h*mean([max(abs(dh_cn_exp)),max(abs(dh_cm_exp))]));  % AoA amplitude [deg]

%% Find equivalent AoA
for coef=1:2
    if coef == 1
        dh = dh_cn_exp;
    elseif coef == 2
        dh = dh_cm_exp;
    end
    % Adjust first and last entries of plunge array
    dh(1) = 0; dh(end) = 0;
    % Nondimensional plunge amplitude []
    A = max(abs(dh));
    % Cycle angle [deg]
    time = asind(dh/A);
    % Adjust cycle angle
    [~,ind_90] = max(dh);
    [~,ind_270] = min(dh);
    time(ind_90+1:ind_270) = 180-time(ind_90+1:ind_270);
    time(ind_270+1:end) = 360+time(ind_270+1:end);
    % Dimensional time
    time_d = T*time/360;
    % Plunge velocity
    hdot = A*2*pi*f*cos(2*pi*f*time_d);
    % Plunge-induced AoA 
    alpha = a0+atand(-hdot*b/U);
    % Set data
    if coef == 1
        alpha_exp_cn = alpha;
        c = cn_exp;
    elseif coef == 2
        alpha_exp_cm = alpha;
        c = cm_exp;
    end
    % Plot
    figure;plot(alpha,c,'k-o');grid; figure;plot(time,alpha,'k-o');grid; 
    % figure;plot(alpha,c,'k-o',alpha_Tyler,cn_Tyler,'b-d');grid
end

%% Save data
new_filename_cn = erase(filename_cn,'_raw');
new_filename_cm = erase(filename_cm,'_raw');
save(new_filename_cn,'cn_exp','alpha_exp_cn');
save(new_filename_cm,'cm_exp','alpha_exp_cm');