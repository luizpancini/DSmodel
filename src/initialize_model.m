function [data,params,x0,y0,RKF_options] = initialize_model(case_now,INPUTS,varargin)

%% Handle inputs
% Default values
n_discard = 3;              % Number of cycles to discard due to initial transients of the model
dt_max = 5e-4;              % Default maximum time step (notice that the NMRS errors are slightly dependent on the time step)
initialize = true;          % Option to initialize states by running one cycle
display_progress = false;   % Option to display algorithm progress
% Input values
if ~isempty(varargin)
    n_discard = varargin{1};
    if length(varargin) > 1
        dt_max = varargin{2};
    end
    if length(varargin) > 2
        initialize = varargin{3};
    end
end

%% Read reference data and simulation parameters from input
[data,params] = read_data(INPUTS,case_now,n_discard);                      

%% Initialize to get consistent set of initial states 
% Set RKK45 ODE solver options
RKF_options.dt_min = 1e-4;                          % Minimum timestep
RKF_options.dt_max = dt_max;                        % Maximum timestep
RKF_options.display_progress = display_progress;    % Option to display progress
RKF_options.RKFtol = 1e-7;                         % Default absolute tolerance for difference between 4th and 5th order RKF approximations
RKF_options.RKF_it_max = 10;                        % Default maximum number of iterations for RKF approximations
% Solve ODEs for initial time
switch INPUTS.model
    case "BL"
        % Set model's number of states
        params.N = 14;    
        % Initialize states and outputs arrays
        x0 = zeros(params.N,1); x0(10:12) = 1; y0 = nan(48,1);
        % Get initial states
        if initialize
            [~,~,x,y] = BL_RKF45([0; params.t_cycle],x0,y0,params,RKF_options);
            x0 = x(:,end); y0 = y(:,end);
        end
    case "BLT"
        % Set model's number of states
        params.N = 8;    
        % Initialize states and outputs arrays
        x0 = zeros(params.N,1); x0(4:6) = 1; y0 = nan(36,1);
        % Get initial states
        if initialize
            [~,~,x] = BLT_RKF45([0; params.t_cycle],x0,y0,params,RKF_options);
            x0 = x(:,end);
        end   
    case {"BLS","BLSLR"}
        % Set model's number of states
        params.N = 9;   
        % Initialize states and outputs arrays
        x0 = zeros(params.N,1); x0(8:9) = 1; y0 = nan(20,1);
        % Set linearized reattachment option for Sheng's model
        if INPUTS.model == "BLSLR", params.lin_reat = true; else, params.lin_reat = false; end
        % Get initial states
        if initialize
            [~,~,x] = BLS_RKF45([0; params.t_cycle],x0,y0,params,RKF_options);
            x0 = x(:,end);
        end
    case "BLO"
        % Set model's number of states
        params.N = 12;   
        % Initialize states and outputs arrays
        x0 = zeros(params.N,1); x0([10,12]) = 1; y0 = nan(24,1);
        % Get initial states
        if initialize
            [~,~,x] = BLO_RKF45([0; params.t_cycle],x0,y0,params,RKF_options);
            x0 = x(:,end);
        end
    case "BLG"
        % Set model's number of states
        params.N = 11;  
        % Initialize states and outputs arrays
        x0 = zeros(params.N,1); x0(10) = 1; y0 = nan(24,1);
        % Get initial states
        if initialize
            [~,~,x] = BLG_RKF45([0; params.t_cycle],x0,y0,params,RKF_options);
            x0 = x(:,end);
        end
    case "BLgust"
        % Set model's number of states
        params.N = 17+params.N_g;    
        % Initialize states and outputs arrays
        x0 = zeros(params.N,1); x0(10:12) = 1; y0 = nan(36,1);
        % Modification for breakpoint of separation angle
        params.alpha1_0 = 17*pi/180;
        % Get initial states
        if initialize && params.t_init > 0
            [~,~,x] = BLgust_RKF45([-params.t_init; 0],x0,y0,params,RKF_options);
            x0 = x(:,end);
        end    
    case "BLTgust"
        % Set model's number of states
        params.N = 9+params.N_g;
        % Initialize states and outputs arrays
        x0 = zeros(params.N,1); x0(4:6) = 1; y0 = nan(37,1);
        % Modification for breakpoint of separation angle
        params.alpha1_0 = 17*pi/180;
        % Get initial states
        if initialize && params.t_init > 0
            [~,~,x] = BLTgust_RKF45([-params.t_init; 0],x0,y0,params,RKF_options);
            x0 = x(:,end);
        end
    case "BLOgust"
        % Set model's number of states
        params.N = 12+params.N_g;   
        % Initialize states and outputs arrays
        x0 = zeros(params.N,1); x0([10,12]) = 1; y0 = nan(25,1);
        params.alpha1_0 = 17*pi/180; params.c_n1 = 1.66; 
        % Get initial states
        if initialize && params.t_init > 0
            [~,~,x] = BLOgust_RKF45([-params.t_init; 0],x0,y0,params,RKF_options);
            x0 = x(:,end);
        end  
    case "BLTflap"
        % Set model's number of states
        params.N = 10;    
        % Initialize states and outputs arrays
        x0 = zeros(params.N,1); x0(4:6) = 1; y0 = nan(39,1);
        % Get initial states
        if initialize
            [~,~,x] = BLTflap_RKF45([0; params.t_cycle],x0,y0,params,RKF_options);
            x0 = x(:,end);
        end
    case "BLtvfs"
        % Set model's number of states
        params.N = 17;  
        % Initialize states and outputs arrays
        x0 = zeros(params.N,1); x0(2:4) = 1; y0 = nan(37,1);
%         % For recurrence algorithm:
%         params.N = 6; y0 = nan(50,1);
%         % For convolution algorithm:
%         params.N = 6; 
    case "BLTtvfs"
        RKF_options.RKFtol = 1e-7;
        % Set model's number of states
        params.N = 8; 
        % Initialize states and outputs arrays
        x0 = zeros(params.N,1); x0(2:4) = 1; y0 = nan(39,1);  
    case "BLhargen"
        % Set model's number of states
        params.N = 22+params.N_g;    
        % Initialize states and outputs arrays
        x0 = zeros(params.N,1); x0(2:4) = 1; y0 = nan(62,1);
        % Increase RKF tolerance for cases of time-varying freestream
        if any(params.U_1), RKF_options.RKFtol = 1e-7; end
        % Check initialization option
        if isfield(data,'source') && data.source == "gust_test" && params.t_init == 0
            initialize = false;
        end
        % Get initial states
        if initialize
            if params.t_init > 0
                tspan = [-params.t_init; 0]; %params.alpha1_0n = 17.0*pi/180;
            else
                tspan = [0; params.t_cycle];
            end
            [~,~,x] = BLhargen_RKF45(tspan,x0,y0,params,RKF_options);
            x0 = x(:,end);
        end 
   case "BLThargen"
        % Set model's number of states
        params.N = 11+params.N_g;    
        % Initialize states and outputs arrays
        x0 = zeros(params.N,1); x0(2:4) = 1; y0 = nan(62,1);
        % Increase RKF tolerance for cases of time-varying freestream
        if any(params.U_1), RKF_options.RKFtol = 1e-7; end
        % Check initialization option
        if isfield(data,'source') && data.source == "gust_test" && params.t_init == 0
            initialize = false;
        end
        % Get initial states
        if initialize
            if params.t_init > 0
                tspan = [-params.t_init; 0]; %params.alpha1_0 = 17.5*pi/180;
            else
                tspan = [0; params.t_cycle];
            end
            [~,~,x] = BLThargen_RKF45(tspan,x0,y0,params,RKF_options);
            x0 = x(:,end);
        end 
end

end