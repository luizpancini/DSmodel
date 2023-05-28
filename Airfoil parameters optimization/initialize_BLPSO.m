function [Aim,bestAimIter,bestGlobalAim,bestOwnAim,bestOwnPos,bestInformantsAim,bestInformantsPos,bestMemInformantsAim,bestMemInformantsPos,bestExpInformantsAim,bestExpInformantsPos,informants,cumStallIt,FE,moveMethodFlag,N_c0,N_i0,N_max,N_min,N_e_range,N_m_range,noiseFlag,P,V,topologyFlag,axes1,axes_size,lw,ms] = initialize_BLPSO(D,EM_swarms,initializePos,initializeVel,LB,maxFE,moveMethod,N,N_c,N_i,N_e,N_m,num_workers_parfor,topology,UB,Vmax,run_parfor,set_GP,specified_filename)

% Arrays and TFs
Aim = nan(N,1);                         % Array for particles' aims
bestAimIter = nan(round(maxFE/N),2);    % Array for iteration's best aim (1st column) and which particle it belongs to (2nd column)
bestGlobalAim = 1e9;                    % Best aim ever among all particles
bestOwnAim = 1e9*ones(N,1);             % Particles' current best personal aim
bestOwnPos = zeros(N,D);                % Particles' best personal position
bestInformantsAim = zeros(N,1);         % Array for best informants' aims in the single swarm case
bestInformantsPos = zeros(N,D);         % Array for best informants' positions in the single swarm case
bestMemInformantsPos = zeros(N_m,D);    % Array for positions of the best informants of the memory swarm 
bestMemInformantsAim = zeros(N_m,1);    % Array for aims of the best informants of the memory swarm 
bestExpInformantsPos = zeros(N_e,D);    % Array for positions of the best informants of the explorer swarm 
bestExpInformantsAim = zeros(N_e,1);    % Array for aims of the best informants of the explorer swarm 
informants = [];                        % Informants array
FE = 0;                                 % Number of function evaluations

% Variables for explorer/memory swarms
if EM_swarms
    N_m_range = 1:N_m;                  % Range of particles in the memory swarm
    N_e_range = N_m+1:N;                % Range of particles in the explorer swarm
else
    N_m_range = [];
    N_e_range = [];
end

% Adaptation variables
N_max = min(200,3*N);                   % Maximum variable swarm size (thrice the initial value, limited to 200)
N_min = max(D+1,N/2);                   % Minimum variable swarm size (half the initial value, limited to D+1)
N_i0 = N_i;                             % Initial / limiting number of informants 
N_c0 = N_c;                             % Initial / limiting number of communication channels
cumStallIt = 0;                         % Cummulative counter of stalled iterations for inertia weight adaptation 

% Initial positions
switch initializePos
    case "random"
        P = LB + rand(N,D) .* (UB-LB);
    case "hemmersly"
        P = hemmersly(N,D,LB,UB);
    otherwise
        error("'" + initializePos + "' position initialization method unavailable");
end

% Initial velocities
switch initializeVel
    case "null"
        V = zeros(N,D);
    case "rand_lim"
        V = -Vmax + 2*rand(N,D) .* Vmax;
    case "rand_shift"
        V = (LB-P) + rand(N,D) .* (UB-LB);
    otherwise
        error("'" + initializeVel + "' velocity initialization method unavailable");
end

% Load guess particle's parameters
if exist(specified_filename,'file') && set_GP > 0
    L = load(specified_filename,'bestGlobalPos','P','V');
    if set_GP == 1
        GP = L.bestGlobalPos;
        GV = [];
    else
        if isfield(L,'P'), GP = L.P; else, GP = L.bestGlobalPos; end
        if isfield(L,'V'), GV = L.V; else, GV = []; end
        if size(GP,2) ~= D
            error("Number of optimization variables from loaded particles is different than the current one")
        end
    end
end

% Set guess particles' positions and velocities
if set_GP > 0
    N_loaded = size(GP,1);
    Np = min([N_loaded, N]);
    P(1:Np,:) = GP;
    if ~isempty(GV)
        V(1:Np,:) = GV;
    end
end

% Flags for movement method 
switch moveMethod
    case "uniform"
        moveMethodFlag = 1;
    case {"pivots","noisy_pivots"}
        moveMethodFlag = 2;
end

% Set noise flag for movement method
if contains(moveMethod,"noisy")
    noiseFlag = 1;
else
    noiseFlag = 0;
end

% Flags for topology method
switch topology
    case "fixed_N_randi"
        topologyFlag = 1;
    case "adaptive_N_randi"
        topologyFlag = 2;    
    case "ring"
        topologyFlag = 3;
    case "SL_ring"
        topologyFlag = 4; 
    case "VR_ring"
        topologyFlag = 5; 
end

% Initialize plot
set(0,'DefaultTextInterpreter','tex')
set(0,'DefaultLegendInterpreter','tex')
axes_size = 20; lw = 0.75; ms = 3;
figure1 = figure('InvertHardcopy','off','Color',[1 1 1]); figure1.PaperSize = [6 4.5]; figure1.Name = "Aim vs iterations";
axes1 = axes('Parent',figure1,'FontSize',axes_size,'FontName','times new roman');
hold(axes1,'on'); grid on; axes1.YScale = 'log';
xlabel('Iteration','FontWeight','normal','FontSize',axes_size);
ylabel('Aim','FontWeight','normal','FontSize',axes_size);

% Initialize parallel pool
if run_parfor
    p = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(p)
        parpool('local',num_workers_parfor);
    end
else
    delete(gcp('nocreate'))
end
    
end