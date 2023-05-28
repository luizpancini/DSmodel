clc
clear 
close all
addpath('../functions')

%% Dynamic stall model inputs
% Select source of reference data
BL_INPUTS.source = "OSU";
% Cases
% BL_INPUTS.cases = [372 373 374 384 385 386 396 397 398 408 409 410 420 421 422 432 433 434]; % M = 0.075     
% BL_INPUTS.cases = [375 376 377 387 388 389 399 400 401 411 412 413 423 424 425 435 436 437]; % M = 0.100 
% BL_INPUTS.cases = [378 379 380 390 391 392 402 403 404 414 415 416 426 427 428 438 439 440]; % M = 0.125
% BL_INPUTS.cases = [381 382 383 393 394 395 405 406 407 417 418 419 429 430 431 441 442 443]; % M = 0.145
BL_INPUTS.cases = 372:443; % All cases
% Select model 
BL_INPUTS.model = "BL"; BL_INPUTS.gust_ind = "";
% Airfoil from cases
BL_INPUTS.airfoil = 'S809';
% Normalized weights for each of the specified cases to be optimized (sum of weights must be equal to 1)
% cases_weights = ones(length(BL_INPUTS.cases),1)/length(BL_INPUTS.cases);
% BL_INPUTS.cases_weights = repmat([1 1 2],1,6)/length(BL_INPUTS.cases);
BL_INPUTS.cases_weights = repmat([1 1 2],1,24)/length(BL_INPUTS.cases);
% Weights for aerodynamic coefficients optimization
BL_INPUTS.coefs_weights.cn = 2;
BL_INPUTS.coefs_weights.cm = 2;
BL_INPUTS.coefs_weights.cc = 2;
% Weights for dynamic stall events
BL_INPUTS.events_weights.ds = 0.5;    % Stall onset angle (angle of maximum c_c)
BL_INPUTS.events_weights.re = 1.0;    % Reattachment (angle and c_n value)
BL_INPUTS.events_weights.le = 0.25;   % Loads extrema (angle and coefficients values)

%% Set list of configurations (PSO algorithm options) to try
configurations_list = [23];

%% Loop over runs
% Number of runs on each configuration
N_runs = 1;
% Loop over runs
for run=1:N_runs
    % Loop over configurations
    for config=configurations_list(:).'
        % Get options from defined configuration
        PSO_options = BLPSO_configurations(BL_INPUTS.model,config);
        % PSO iterations
        [bestGlobalAim,bestGlobalPos,P,V] = BLPSO_speedcore(BL_INPUTS,PSO_options);
        % Plot results for overall best
        fprintf('\n Best overall aim was %.4e \n', bestGlobalAim);
        [~,bestParams,E] = optimize_BL(bestGlobalPos,BL_INPUTS,1);
        % Save best on current run
        next_filename = nextname('gBest','_001','.mat');
        save(next_filename,'bestGlobalPos','bestGlobalAim','bestParams','BL_INPUTS','PSO_options','E','config','P','V');
    end
end
