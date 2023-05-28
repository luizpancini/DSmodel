clc
clear
close all
addpath('../functions')

%% Dynamic stall model inputs
% Select source 
BL_INPUTS.source = "gust_test";
% Cases
BL_INPUTS.cases = [14,15,16,17,19,20]; 
% Select model 
BL_INPUTS.model = "BLgust"; 
BL_INPUTS.gust_ind = "BR-C";
% Airfoil from cases
BL_INPUTS.airfoil = 'NACA0006';

%%
% Option to set initial particles in the swarm (Set as 0 to not load a guess particle,
% set as 1 to load the best from a specified file, or 2 to load all from a
% specified file. If no such file exists, do nothing)
set_GP = 0;
specified_filename = 'gBestMat_007.mat';
% Number of particles
N = 50; 
% Number of dimensions (optimization variables)
D = 18; 
% Maximum number of iterations
MaxIt = 150;
% Maximum number of iteration without improvement to detect as stagnation
MaxStallIt = 30;
% Bounds according to airfoil
[LB,UB] = set_BLPSO_cgust_bounds(N,BL_INPUTS.airfoil);
% Define optimization function
fun = @(P) optimize_BL_cgust(P,BL_INPUTS);

%%
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
options.Display = 'iter';
% options.HybridFcn = @fmincon;
options.MaxIterations = MaxIt;
options.PlotFcn = @pswplotbestf;
% options.SwarmSize = N;
options.UseParallel = true;
options.MaxStallIterations = MaxStallIt;
if set_GP > 0
    options.InitialSwarmMatrix = GP;
end
% PSO
[bestGlobalPos,bestGlobalAim,exitflag,output] = particleswarm(fun,D,LB,UB,options);

%% Save and plot
% Pack algorithm data into PSO structure
PSO = variables2struct(struct,exitflag,LB,MaxIt,N,options,output,D,set_GP,UB);
% Save best on current run
next_filename = nextname('gBestMat','_001','.mat');
save(next_filename,'bestGlobalPos','bestGlobalAim','BL_INPUTS','PSO');
% Plot results for overall best
optimize_BL_cgust(bestGlobalPos,BL_INPUTS,1);
