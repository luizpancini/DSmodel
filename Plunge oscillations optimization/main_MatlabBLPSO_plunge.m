clc
clear
% close all
addpath('../functions')

%% Inputs
% Select source, airfoil and model
source = "other";
airfoil = 'Vertol-23010';
model = "BL"; 
% Cases at each airspeed available
M_02_cases = [];                      % M = 0.2 
M_04_cases = 13:20;                   % M = 0.4 
M_06_cases = [];                      % M = 0.6
% Number of runs for each airspeed
N_runs_M_02 = 0;
N_runs_M_04 = 2;
N_runs_M_06 = 0;
% Optimization weights
cases_weights = [1, 1, 3, 3, 3, 3, 3, 3];
% cases_weights = [1, 3, 3, 3];
coefs_weights = [1, 1, 1];            % c_n, c_m, c_c
%%% Algorithm options
% Option to set initial particles in the swarm (Set as 0 to not load a guess particle, set as 1 to load the best from a specified file, or 2 to load all from a specified file. If no such file exists, do nothing)
set_GP = 1; specified_filename = 'gBestMat_005.mat';
% Number of particles
N = 60; 
% Maximum number of iterations
MaxIt = 150;
% TF to use parallel
UseParallel = false;

%% Setup inputs
% Model
BL_INPUTS.source = source;
BL_INPUTS.airfoil = airfoil;
BL_INPUTS.model = model; 
BL_INPUTS.gust_ind = "";
% Weights for aerodynamic coefficients 
BL_INPUTS.coefs_weights.cn = coefs_weights(1);
BL_INPUTS.coefs_weights.cm = coefs_weights(2);
BL_INPUTS.coefs_weights.cc = coefs_weights(3);
% Cases and number of runs definition
N_cases = [N_runs_M_02,N_runs_M_04,N_runs_M_06];
N_cases_cumsum = cumsum(N_cases);
N_runs = sum(N_cases);
% Optimization variables' bounds according to airfoil
[D,LB,UB] = set_BLPSO_plunge_bounds(N,BL_INPUTS.airfoil);

%% PSO setup 
% Load guess particle's parameters
if exist(specified_filename,'file') && set_GP > 0
    load(specified_filename,'bestGlobalPos','PSO');
    switch set_GP
        case 1
            GP = bestGlobalPos;
        case 2
            GP = PSO.P;
            clear PSO;
    end
    if size(GP,2) ~= D
        error("Number of optimization variables from loaded particles is different than the current one")
    end
end
% Optimization options
options = optimoptions('particleswarm');
if set_GP > 0, options.InitialSwarmMatrix = GP; end
options.Display = 'iter';
options.MaxIterations = MaxIt;
options.PlotFcn = @pswplotbestf;
options.SwarmSize = N;
options.UseParallel = UseParallel;
% options.HybridFcn = @fmincon;

%% Loop over runs
for run=1:N_runs
    % Current cases
    f = find(run>N_cases_cumsum,1,'last');
    if isempty(f)
        cases_now = M_02_cases;
    else
        switch f
            case 1
                cases_now = M_04_cases;
            case 2
                cases_now = M_06_cases;
        end
    end
    % Cases to run now and normalized weights
    BL_INPUTS.cases = cases_now;
    BL_INPUTS.cases_weights = cases_weights/sum(cases_weights); 
    % Define optimization function
    fun = @(P) optimize_BL_plunge(P,BL_INPUTS);
    % PSO
    [bestGlobalPos,bestGlobalAim,exitflag,output] = particleswarm(fun,D,LB,UB,options);
    % Pack algorithm data into PSO structure
    PSO = variables2struct(struct,exitflag,LB,MaxIt,N,options,output,D,set_GP,UB);
    % Plot results for overall best
    [E_sum,E,bestParams] = optimize_BL_plunge(bestGlobalPos,BL_INPUTS,1);
    % Save best on current run
    next_filename = nextname('gBestMat','_001','.mat');
    save(next_filename,'bestGlobalPos','bestGlobalAim','BL_INPUTS','PSO','E','bestParams');
end