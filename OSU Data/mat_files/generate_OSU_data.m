clc
clear

%% Set file list
file_list = {'C5l75_s809_8'    'C5m75_s809_8'    'C5h75_s809_8'
             'C5l100_s809_8'   'C5m100_s809_8'   'C5h100_s809_8'
             'C5l125_s809_8'   'C5m125_s809_8'   'C5h125_s809_8'
             'C5l140_s809_8'   'C5m140_s809_8'   'C5h140_s809_8'
             'C5l75_s809_14'   'C5m75_s809_14'   'C5h75_s809_14'
             'C5l100_s809_14'  'C5m100_s809_14'  'C5h100_s809_14'
             'C5l125_s809_14'  'C5m125_s809_14'  'C5h125_s809_14'
             'C5l140_s809_14'  'C5m140_s809_14'  'C5h140_s809_14'
             'C5l75_s809_20'   'C5m75_s809_20'   'C5h75_s809_20'
             'C5l100_s809_20'  'C5m100_s809_20'  'C5h100_s809_20'
             'C5l125_s809_20'  'C5m125_s809_20'  'C5h125_s809_20'
             'C5l140_s809_20'  'C5m140_s809_20'  'C5h140_s809_20'
             'C10l75_s809_8'   'C10m75_s809_8'   'C10h75_s809_8'
             'C10l100_s809_8'  'C10m100_s809_8'  'C10h100_s809_8'
             'C10l125_s809_8'  'C10m125_s809_8'  'C10h125_s809_8'
             'C10l140_s809_8'  'C10m140_s809_8'  'C10h140_s809_8'
             'C10l75_s809_14'  'C10m75_s809_14'  'C10h75_s809_14'
             'C10l100_s809_14' 'C10m100_s809_14' 'C10h100_s809_14'
             'C10l125_s809_14' 'C10m125_s809_14' 'C10h125_s809_14'
             'C10l140_s809_14' 'C10m140_s809_14' 'C10h140_s809_14'
             'C10l75_s809_20'  'C10m75_s809_20'  'C10h75_s809_20'
             'C10l100_s809_20' 'C10m100_s809_20' 'C10h100_s809_20'
             'C10l125_s809_20' 'C10m125_s809_20' 'C10h125_s809_20'
             'C10l140_s809_20' 'C10m140_s809_20' 'C10h140_s809_20'};
file_list = {'C5m75_s809_8'};

%% Constants
authors{1} = 'NREL/OSU';
airfoil = "S809";
b = 0.457/2; % Semi-chord
a_inf = 343; % Sound speed [m/s]

% Interpolation method
interp_method = 'pchip'; % Choose 'spline' or 'pchip'

% Plotting option (TF)
plotting = 1;

%% Loop over files         
for f = 1:numel(file_list)
    if plotting, close all; end
    %% Filename
    filename = file_list{f};
    %% Get file path
    filepath = "../S809/" + filename + ".txt";
    text = fileread(filepath);
    %% Read test conditions data
    fid = fopen(filepath);
    RUN = str2double(regexpi(text,'(?<=RUN  )\d*','match'));
    U = 0.3051*str2double(regexpi(text,'(?<=TUNNEL AIRSPEED  =  )\d+\.?\d*','match'));
    if isempty(U)
        U = 0.3051*str2double(regexpi(text,'(?<=TUNNEL AIRSPEED  = )\d+\.?\d*','match'));
    end
    Re = 1e6*str2double(regexpi(text,'(?<=REYNOLDS NUMBER  =  )\d+\.?\d*','match'));
    omega_Hz = str2double(regexpi(text,'(?<=OSCILLATOR FREQUENCY =  )\d+\.?\d*','match'));
    M = U/a_inf;
    k = 2*pi*omega_Hz*b/U;
    T = 1./omega_Hz;
    beta = sqrt(1-M^2);
    fclose(fid);
    %% Get coefficients data
    OSU_data = importdata(filepath);
    time_s = OSU_data.data(:,2);                % Time (s)
    alpha = OSU_data.data(:,3);                 % AoA (deg)
    cl = OSU_data.data(:,4);
    cd = OSU_data.data(:,5);
    cm = OSU_data.data(:,6);
    a_0 = (max(alpha)+min(alpha))/2*pi/180;
    a_1 = (max(alpha)-min(alpha))/2*pi/180;
    cn = cl.*cosd(alpha)+cd.*sind(alpha);
    cc = cl.*sind(alpha)-cd.*cosd(alpha);
    %% Find indices where new cycles begin
    ind = 1; n = 2;
    for i=1:length(time_s)
        if time_s(i)/T > n-1
            ind(n) = i;
            n = n+1;
        end
    end
    n = n-1; % Total number of cycles (including partial ones)     
    %% Cycle angle [deg] at which the data begins
    if alpha(2) >= alpha(1) % Starting on upstroke
        time_init = 180/pi*real(asin((alpha(1)*pi/180-a_0)/a_1));
    else                    % Starting on downstroke
        time_init = 180-180/pi*real(asin((alpha(1)*pi/180-a_0)/a_1));
    end
    %% Transform time to cycle angle [deg]
    time = time_init + 360*rem(time_s/T,1);
    if time(1) < -270, time = time+360; end  % Set initial time greater than -270 deg
    if time(end) > 270, time = time-360; end % Set final time smaller than 270 deg
    %% Get initial interpolated coefficients and AoA separated by cycle
    time_ref = linspace(time(1),time(1)+360,length(time)+1);
    alpha_ref = 180/pi* (a_0 + a_1*sind(time_ref));
    if plotting, figure; plot(time_ref,alpha_ref,'k-'); hold on; grid; end
    for i=1:n-1
        cycle_range = ind(i):ind(i+1)-1;
        time_cycles{i} = time(cycle_range);
        alpha_cycles{i} = alpha(cycle_range);
        cl_cycles{i} = cl(cycle_range);
        cd_cycles{i} = cd(cycle_range);
        cm_cycles{i} = cm(cycle_range);
        cn_cycles{i} = cn(cycle_range);
        cc_cycles{i} = cc(cycle_range);
        alpha_interp(:,i) = interp1(time_cycles{i},alpha_cycles{i},time_ref,interp_method);
        if plotting
            plot(time_cycles{i},alpha_cycles{i},'-o');
        end
    end
    % Mean interpolated AoA
    alpha_interp = sum(alpha_interp,2)/(n-1);
    if plotting
        plot(time_ref,alpha_interp,'-');
    end
    %% Find time shift so that mean interpolated AoA and the reference AoA input match at the cycle angle = 0 deg (begin of cycle)
    % Index at which interpolated AoA is closest to mean angle
    [~,ind_interp] = min(abs(alpha_interp*pi/180-a_0));
    % Cycle angle shift between reference and interpolated AoA curves
    t_shift = rem(time_ref(ind_interp),180);  
    % Limit to the range [-90,90]
    if t_shift > 90, t_shift = t_shift-180; end 
    if t_shift < -90, t_shift = -t_shift+180; end
    % True initial time for matching
    time_init = time_init - t_shift; 
    %% Repeat calculations with true initial time
    % Transform time to cycle angle [deg]
    time = time_init + 360*rem(time_s/T,1);
    if time(end) > 270, time = time-360; end
    if time(1) < -270, time = time+360; end
    % Get coefficients interpolated coefficients and AoA separated by cycle
    time_ref = linspace(-90,270,length(time)+1);
    alpha_ref = 180/pi* (a_0 + a_1*sind(time_ref));
    if plotting, figure; plot(time_ref,alpha_ref,'k-');  hold on; grid; xlim([-90 270]); ax = gca; ax.XTick = [-90 0 90 180 270]; end
    for i=1:n-1
        % Find current cycle's range
        cycle_range = ind(i):ind(i+1)-1;
        % Find first index such that cycle angle is greater than -90
        ind_gm90 = find(time(cycle_range)>=-90,1,'first')-1;       
        % Find index such that cycle angle is closest to 270
        [~,ind_270] = min(abs(time(cycle_range)-270));
        % If ind_gm90 is zero and ind_270 is not the last in the cycle, then
        % the cycle angle must be shifted to the left
        if ind_gm90 == 0 && ind_270 ~= length(cycle_range)
            ind_g270 = find(time(cycle_range)>=270,1,'first');
            ind_shift = length(cycle_range)-ind_g270+1;
            shifted_range = 1:ind_shift;
            delta_angle = -360;
        % Otherwise, it must be shifted to the right    
        else
            ind_shift = -ind_gm90;
            shifted_range = length(cycle_range)-ind_gm90+1:length(cycle_range);
            delta_angle = 360;
        end
        % Shifted values for each cycle
        time_cycles{i} = circshift(time(cycle_range),ind_shift); 
        time_cycles{i}(shifted_range) = time_cycles{i}(shifted_range)+delta_angle;
        alpha_cycles{i} = circshift(alpha(cycle_range),ind_shift); 
        cl_cycles{i} = circshift(cl(cycle_range),ind_shift); 
        cd_cycles{i} = circshift(cd(cycle_range),ind_shift); 
        cm_cycles{i} = circshift(cm(cycle_range),ind_shift); 
        cn_cycles{i} = circshift(cn(cycle_range),ind_shift); 
        cc_cycles{i} = circshift(cc(cycle_range),ind_shift); 
        % Set additional last value equal to the first for each cycle,
        % so that interpolated values do not differ much at begin and end
        time_cycles{i}(end+1) = time_cycles{i}(1)+360;
        alpha_cycles{i}(end+1) = alpha_cycles{i}(1);
        cl_cycles{i}(end+1) = cl_cycles{i}(1);
        cd_cycles{i}(end+1) = cd_cycles{i}(1);
        cm_cycles{i}(end+1) = cm_cycles{i}(1);
        cn_cycles{i}(end+1) = cn_cycles{i}(1);
        cc_cycles{i}(end+1) = cc_cycles{i}(1);
        % Plots
        if plotting
            alpha_interp(:,i) = interp1(time_cycles{i},alpha_cycles{i},time_ref,interp_method);
            plot(time_cycles{i},alpha_cycles{i},'-o');
%             cl_interp(:,i) = interp1(time_cycles{i},cl_cycles{i},time_ref,'spline');
%             plot(time_cycles{i},cl_cycles{i},'o');
        end
    end
    % Mean interpolated AoA
    if plotting
        alpha_interp = sum(alpha_interp,2)/(n-1);
        plot(time_ref,alpha_interp,'-');
    end
    %% Save data
    filepath_save = sprintf('OSU_%d.mat', RUN);
    save(filepath_save,'authors','airfoil','b','a_inf','Re','M','U','beta','k','a_0','a_1','time','alpha','cl','cd','cm','cn','cc','time_cycles','alpha_cycles','cl_cycles','cd_cycles','cm_cycles','cn_cycles','cc_cycles');
end
