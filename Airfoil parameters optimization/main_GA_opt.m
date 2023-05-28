clc
clear
close all
addpath('../functions')

%%
% Current optimization parameters
% #1: alpha_ss
% #2: alpha_ds0
% #3: r0
% #4: Ta

% Set lower bounds for parameters:
lb = [14.8; 18.0; 0.001; 2.0];

% Set upper bounds for parameters:
ub = [15.3; 21.0; 0.006; 5.0];

% Start the optimization
options = optimoptions('ga','PlotFcn','gaplotbestf','PopulationSize',25,'UseParallel',1);
[p_best,best_E,exitflag,output,population,scores] = ga(@optimize_ds,4,[],[],[],[],lb,ub,[],options);