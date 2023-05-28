% clc
clear
% close all
addpath('../functions')

%% Inputs
% Select source, airfoil and model
source = "GU";
airfoil = 'NACA0018';
model = "BL"; 
%%%%%%%%%%%%%%%%%%% Cases to run for NACA0012 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cases = [7019,7021,7023,7101,7104,7108,7110:7114,7117:7121,7200,7202,7205,7207,7212,7214,7216,7222,7300,7302,7305,9210,9213,9214,9221:9223,9302,9307,10022,10104,10105,10108,10113,10114,10117,10120,10123,10202,10203,10204,10207,10208,10211,10212,10218,10221,10222,10303,10305,10309];  % M = 0.3  
% cases = [12109,9208,9218]; % M = 0.28
% cases = [9203,13303]; % M = 0.25
% cases = [9202,13308]; % M = 0.215
% cases = [8220,8222,8306,9022,9101,9106,9110,9112,9118]; % M = 0.185
% cases = 11014311:10:11014461; % M = 0.155
% cases = 11000041; % M = 0.155 (STATIC)
% cases = [11000011,11011962:10:11012032,11012052:10:11012282,11012302:10:11012442,11012622:10:11012862,11012891:10:11012961,11012981:10:11013191,11013211:10:11013351,11013371:10:11013471,11013491:10:11013741,11013761:10:11014121]; % M = 0.117
% cases = [11011962:10:11012032,11012052:10:11012282,11012302:10:11012442,11012622:10:11012862];
% cases = [11012062,11012102,11012152,11012192,11012242,11012282,11012362,11012412,11012692];
% cases = [11012931:10:11013041,11013091:10:11013191,11013251:10:11013351,11013411:10:11013471,11013491:10:11013511,11013561:10:11013661,11013711:10:11013821,11013881:10:11013971,11014031:10:11014121];
% cases = [11013091:10:11013191,11013261:10:11013351,11013421:10:11013471,11013491:10:11013511,11013581:10:11013661,11013741,11013761:10:11013821,11013901:10:11013971,11014051:10:11014121];
% cases = [8210,8214,13021];    % M = 0.11
% cases = [11014141:10:11014291];   % M = 0.078
% cases = 11000031;   % M = 0.078 (STATIC)
% cases = [8114,8116,8118,8123,8203]; % M = 0.07
% cases = [8019,8021,8023,8102,8104,8106]; % M = 0.035
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Cases to run for NACA0018 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cases = [06019201,06019221,06013741:10:06013791,06013811:10:06013851]; % M = 0.150
cases = [06012631:10:06012641,06012691:10:06012781,06012851:10:06012921,06012991:10:06013061,06013131:10:06013201,06013271:10:06013341,06013411:10:06013481,06014201:10:06014211,06014251:10:06014261,06019081,06019101:20:06019181,06019381:20:06019421]; % M = 0.12 
% cases = [06013631:10:06013731]; % M = 0.080
% cases = [06011001:20:06011081,06019001,06019021]; % M = 0.062
% Number of runs 
N_runs = 1;
% Optimization weights            
coefs_weights = [1, 1, 1];          % c_n, c_m, c_c
% Error calculation
error_standard = 'RMSE';    % RMSE or NRMSE
%%% Algorithm options
% Option to set initial particles in the swarm (Set as 0 to not load a guess particle, set as 1 to load the best from a specified file, or 2 to load all from a specified file. If no such file exists, do nothing)
set_GP = 1; specified_filename = 'gBestMat_044.mat';
% Number of particles
N = 90; 
% Maximum number of iterations
MaxIt = 100;
% TF to use parallel
UseParallel = 1;

%% Process inputs
% Model
BL_INPUTS.source = source;
BL_INPUTS.airfoil = airfoil;
BL_INPUTS.model = model; 
BL_INPUTS.gust_ind = 'K';
% Cases to run 
BL_INPUTS.cases = cases;
% Weights for aerodynamic coefficients and error standard
BL_INPUTS.coefs_weights = coefs_weights;
BL_INPUTS.error_standard = error_standard;
% Optimization variables' bounds according to airfoil
[D,LB,UB] = set_BLPSO_bounds(N,BL_INPUTS.airfoil);
% Define optimization function
fun = @(P) optimize_BL(P,BL_INPUTS);
    
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
    % PSO
    [bestGlobalPos,bestGlobalAim,exitflag,output] = particleswarm(fun,D,LB,UB,options);
    % Pack algorithm data into PSO structure
    PSO = variables2struct(struct,exitflag,LB,MaxIt,N,options,output,D,set_GP,UB);
    % Plot results for overall best
    [E_sum,bestParams,E] = optimize_BL(bestGlobalPos,BL_INPUTS,1);
    % Save best on current run
    next_filename = nextname('gBestMat','_001','.mat');
    save(next_filename,'bestGlobalPos','bestGlobalAim','BL_INPUTS','PSO','E','bestParams');
end