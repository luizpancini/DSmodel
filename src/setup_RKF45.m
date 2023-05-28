function [display_progress,dtf_steps,c5,c45,dt,dt_min,dt_max,i,tc,ti,tf,RKF_it_max,RKFtol,tp,xp,yp,xdotp,x_i,xdot_i,xdot] = setup_RKF45(tspan,x0,y0,RKF_options)

% Unpack RKF45 options
dt_min = RKF_options.dt_min;
dt_max = RKF_options.dt_max;
RKFtol = RKF_options.RKFtol;
RKF_it_max = RKF_options.RKF_it_max;
display_progress = RKF_options.display_progress;

% RKF45 integration variables
dtf_steps = [0, 1/4, 3/8, 12/13, 1, 1/2];           % Time step fractions
c4 = [25/216 0 1408/2565 2197/4104 -1/5 0]';        % Coefficients for 4th order approximation
c5 = [16/135 0 6656/12825 28561/56430 -9/50 2/55]'; % Coefficients for 5th order approximation
c45 = c5-c4;

% Time variables
ti = tspan(1);              % Initial time
tf = tspan(end);            % Final time
tc = ti;                    % Current time
i = 1;                      % Initialize time step counter
dt = dt_max;                % Initialize time step with maximum allowable

% Pre-allocate
N = length(x0);             % Number of states
s0 = round((tf-ti)/dt_min); % Initial size of the arrays
tp = zeros(s0,1);           % Storage time array
xp = zeros(N,s0);           % Storage states array
yp = zeros(length(y0),s0);  % Storage outputs array
xdotp = zeros(N,s0);        % Storage states' rates array
xdot = zeros(N,6);          % RKF steps' states rates array

% Set initial conditions on current states and states' rates arrays
x_i = x0; 
xdot_i = xdotp(:,1);

% Initialize storage arrays
tp(1) = ti;  
xp(:,1) = x0; 
yp(:,1) = y0;

end