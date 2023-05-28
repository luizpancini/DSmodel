clc
clear
close all

% File to load
filename = 'gBest_022.mat';
load(filename,'bestGlobalPos','PSO','BL_INPUTS');
% BL_INPUTS.cases = 386;
% PSO.cases_weights = 1;

% Bounds according to airfoil
[LB,UB] = set_PSO_bounds(1,'S809');

% Set optimization function and options
ff = @(P) optimize_BL(P,BL_INPUTS,PSO.cases_weights,PSO.coefs_weights,PSO.events_weights);
TolFun = 1e-6; 
MaxFunEvals = 1e3;
options = optimoptions('fmincon','Display','iter','MaxFunEvals',MaxFunEvals,'TolFun',TolFun,'UseParallel',1);

% Find local minimum
[newGBest,fBest,exitflag,output] = fmincon(ff,bestGlobalPos,[],[],[],[],LB,UB,[],options);
