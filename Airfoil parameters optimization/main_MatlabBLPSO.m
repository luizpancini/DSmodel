clc
clear
% close all
addpath('../functions')

%% Inputs
% Select source, airfoil and model
source = "OSU";
airfoil = 'S809';
model = "BL"; 
% Cases at each airspeed available
M_0075_cases = [372 373 374 384 385 386 396 397 398 408 409 410 420 421 422 432 433 434];    % M = 0.075 
M_0100_cases = [375 376 377 387 388 389 399 400 401 411 412 413 423 424 425 435 436 437];    % M = 0.100 
M_0125_cases = [378 379 380 390 391 392 402 403 404 414 415 416 426 427 428 438 439 440];    % M = 0.125
M_0145_cases = [381 382 383 393 394 395 405 406 407 417 418 419 429 430 431 441 442 443];    % M = 0.145
all_cases    = 372:443;                                                                      % All cases 
custom_cases = [429 430 431];
% Number of runs for each airspeed
N_runs_M_0075 = 2;
N_runs_M_0100 = 0;
N_runs_M_0125 = 0;
N_runs_M_0145 = 0;
N_runs_all_cases = 0;
N_runs_custom_cases = 0;
% Optimization weights            
cases_weights = [1, 1, 2];            % More weight on higher reduced frequency cases
coefs_weights = [2, 2, 2];            % c_n, c_m, c_c
events_weights = [0.5, 1.0, 0.25];    % Stall onset, reattachment, loads extrema
%%% Algorithm options
% Option to set initial particles in the swarm (Set as 0 to not load a guess particle, set as 1 to load the best from a specified file, or 2 to load all from a specified file. If no such file exists, do nothing)
set_GP = 0; specified_filename = 'gBestMat_006.mat';
% Number of particles
N = 120; 
% Maximum number of iterations
MaxIt = 120;
% TF to use parallel
UseParallel = true;

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
% Weights for dynamic stall events
BL_INPUTS.events_weights.ds = events_weights(1);    % Stall onset angle (angle of maximum c_c)
BL_INPUTS.events_weights.re = events_weights(2);    % Reattachment (angle and c_n value)
BL_INPUTS.events_weights.le = events_weights(3);    % Loads extrema (angle and coefficient value for c_n and c_m)
% Cases and number of runs definition
N_cases = [N_runs_M_0075,N_runs_M_0100,N_runs_M_0125,N_runs_M_0145,N_runs_all_cases,N_runs_custom_cases];
N_cases_cumsum = cumsum(N_cases);
N_runs = sum(N_cases);
% Optimization variables' bounds according to airfoil
[D,LB,UB] = set_BLPSO_bounds(N,BL_INPUTS.airfoil);

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
        cases_now = M_0075_cases;
    else
        switch f
            case 1
                cases_now = M_0100_cases;
            case 2
                cases_now = M_0125_cases;
            case 3
                cases_now = M_0145_cases;
            case 4
                cases_now = all_cases;
            case 5
                cases_now = custom_cases;
        end
    end
    % Cases to run now and normalized weights
    BL_INPUTS.cases = cases_now;
    BL_INPUTS.cases_weights = repmat(cases_weights,1,length(cases_now)/length(cases_weights))/length(cases_now); 
    % Define optimization function
    fun = @(P) optimize_BL(P,BL_INPUTS);
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