function options = BLPSO_configurations(model,config)

%% Set default options
options = default_options(model);

%% Change options according to configuration
switch config
    case 1  % Default
    case 2  % Test adaptive_N_randi topology
        options.topology = "adaptive_N_randi";
    case 3  % Test ring topology
        options.topology = "ring";
    case 4  % Test SL_ring topology
        options.topology = "SL_ring";
    case 5  % Test VR_ring topology
        options.topology = "VR_ring";
    case 6  % Test pivots moveMethod
        options.moveMethod = "pivots";
    case 7  % Test noisy_pivots moveMethod
        options.moveMethod = "noisy_pivots";
    case 8  % Test adaptiveInertia
        options.adaptiveInertia = 1;
    case 9  % Test adaptiveSwarmSize
        options.adaptiveSwarmSize = 1; 
    case 10 % Test N_i = 4
        options.N_i = 4;
    case 11 % Test EM_swarms
        options.EM_swarms = 1; 
    case 12 % Test EM_swarms with adaptive_N_randi topology
        options.EM_swarms = 1; 
        options.topology = "adaptive_N_randi";
    case 13 % Test EM_swarms with VR_ring topology
        options.EM_swarms = 1; 
        options.topology = "VR_ring";
    case 14 % Test EM_swarms with adaptiveInertia
        options.EM_swarms = 1; 
        options.adaptiveInertia = 1; 
    case 15 % Test EM_swarms with adaptiveInertia and adaptiveSwarmSize
        options.EM_swarms = 1; 
        options.adaptiveInertia = 1;  
        options.adaptiveSwarmSize = 1; 
        options.setNMLOnMaxStalledIt = 1;
    case 16 % Test adaptiveInertia, adaptiveSwarmSize and adaptive_N_randi topology
        options.adaptiveInertia = 1;  
        options.adaptiveSwarmSize = 1; 
        options.topology = "adaptive_N_randi";
    case 17 % Test setNMLOnMaxStalledIt
        options.setNMLOnMaxStalledIt = 1;
    case 18 % Test VR_ring topology, adaptiveSwarmSize, setNMLOnMaxStalledIt and Clerc's coefficients
        options.adaptiveSwarmSize = 1; 
        options.c0 = 1/(2*log(2));
        options.c1 = 1/2 + log(2);
        options.c2 = options.c1;
        options.setNMLOnMaxStalledIt = 1;
        options.topology = "VR_ring";
    case 19 % Test selective informant update with adaptive swarm size
        options.adaptiveSwarmSize = 1;
        options.alwaysUpdateInf = 0;   
    case 20 % Test selective informant update with adaptive swarm size and SL_ring topology
        options.adaptiveSwarmSize = 1;
        options.alwaysUpdateInf = 0;      
        options.topology = "SL_ring";
    case 21 % Test adaptiveInertia, adaptiveSwarmSize and setNMLOnMaxStalledIt
        options.adaptiveInertia = 1;  
        options.adaptiveSwarmSize = 1; 
        options.setNMLOnMaxStalledIt = 1;
%         options.initializePos = "random";
    case 22 % Test adaptiveInertia, adaptiveSwarmSize, setNMLOnMaxStalledIt and SL_ring topology
        options.adaptiveInertia = 1;  
        options.adaptiveSwarmSize = 1; 
        options.setNMLOnMaxStalledIt = 1;
        options.topology = "SL_ring";
    case 23 % Test adaptiveInertia, adaptiveSwarmSize, setNMLOnMaxStalledIt and adaptive_N_randi topology
        options.adaptiveInertia = 1;  
        options.adaptiveSwarmSize = 1; 
        options.setNMLOnMaxStalledIt = 1;
        options.topology = "adaptive_N_randi"; 
    case 24 % Test adaptiveInertia, adaptiveSwarmSize, setNMLOnMaxStalledIt and selective informant update
        options.adaptiveInertia = 1;  
        options.adaptiveSwarmSize = 1; 
        options.setNMLOnMaxStalledIt = 1;
        options.alwaysUpdateInf = 0;
    case 25 % Test adaptiveSwarmSize, setNMLOnMaxStalledIt and selective informant update 
        options.adaptiveSwarmSize = 1; 
        options.setNMLOnMaxStalledIt = 1;
        options.alwaysUpdateInf = 0;    
    case 26
        options.set_GP = 1;                                                 
        options.specified_filename = 'gBest_032.mat'; 
        options.adaptiveInertia = 1;  
        options.adaptiveSwarmSize = 1; 
        options.setNMLOnMaxStalledIt = 1;
    case 27
        options.adaptiveSwarmSize = 1; 
        options.setNMLOnMaxStalledIt = 1;
    case 28
        options.topology = "adaptive_N_randi";
        options.adaptiveInertia = 1;  
        options.N = 100;
        options.setNMLOnMaxStalledIt = 0; 
    case 29
        options.topology = "adaptive_N_randi";
        options.adaptiveInertia = 0;  
        options.N = 100;
        options.setNMLOnMaxStalledIt = 0; 
        options.phi = 2.2;
        options.set_GP = 1;                                                 
        options.specified_filename = 'gBest_050.mat';                       
end

%% Nested functions
    function options = default_options(model)
        switch model
            case "BL"
                options.D = 61;                                             % Number of dimensions (optimization variables)
        end
        % Set defaults
        options.phi = 2.10;                                                 % Constriction coefficient
        options.N = 40+2*round(sqrt(options.D));                            % Swarm size (initial, if adaptive)
        options.adaptiveInertia = 0;                                        % Option for adaptive inertia weight
        options.adaptiveSwarmSize = 0;                                      % Option for adaptive swarm size
        options.alwaysUpdateInf = 1;                                        % Option to always update informants or only if particle's position hasn't improved
        options.c0 = (options.phi-1+sqrt(options.phi^2-2*options.phi))^-1;  % Inertia weight
        options.c1 = options.phi*options.c0;                                % Self-experience weight
        options.c2 = options.phi*options.c0;                                % Social weight
        options.EM_swarms = 0;                                              % Option to set explorer and memory swarms
        options.initializePos = "hemmersly";                                % Initialization method for position. Choose from {"random","hemmersly"}
        options.initializeVel = "rand_shift";                               % Initialization method for velocity. Choose from {"null","rand_lim","rand_shift"}
        options.improveTol = 5e-3;                                          % Tolerance for detection of improvement from one iteration to the next (the new best aim has to be lower than (1-improveTol) times the previous best)
        options.maxFE = 15e3;                                               % Maximum number of function evaluations per run
        options.maxStalledIt = round(options.N/2);                          % Number of iterations without improvement to detect algorithm stall
        options.max_c0 = 1.1;                                               % Upper bound for adaptive inertia weight
        options.min_c0 = 0.1;                                               % Lower bound for adaptive inertia weight
        options.moveMethod = "uniform";                                     % Particle movement method. Choose from {"uniform","pivots","noisy_pivots"}
        options.N_c = ceil(options.N/10);                                   % Number of communication channels, for E/M swarms case (initial, if adaptive)
        options.N_i = 3;                                                    % Number of informants for each particle, including itself, for single swarm case (initial, if adaptive)
        options.N_m = ceil(options.N/2);                                    % Number of particles in the memory swarm (initial, if adaptive)
        options.N_e = options.N-options.N_m;                                % Number of particles in the explorer swarm (initial, if adaptive)
        options.num_workers_parfor = 6;                                     % Number of cores dedicated to parallel run
        options.reboundVel = 1/2;                                           % Rebound velocity ratio for particles going out of bounds (0-1)
        options.run_parfor = 1;                                             % Option to run in parallel mode
        options.set_GP = 0;                                                 % Option to set initial "guess" particles' position and velocity (Set as 0 to not load a guess particle, set as 1 to load the best from a specified file, or 2 to load all from a specified file)
        options.specified_filename = 'gBest_001.mat';                       % Specified filename to load initial "guess" particles
        options.stopOnMaxStalledIt = 0;                                     % Option to stop when maximum number of stalled iterations is reached
        options.setNMLOnMaxStalledIt = 1;                                   % Option to reset a particle into the No Man's Land when maximum number of stalled iterations is reached
        options.topology = "fixed_N_randi";                                 % Topology for particle communication. Choose from {"fixed_N_randi","adaptive_N_randi","ring","SL_ring","VR_ring"}
    end

end