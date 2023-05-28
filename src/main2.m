% Main program
% clc
clear 
% close all

%% Inputs
% Input options:
NASA = [10022];                                                               % NASA experiment #
GU = 11012812;                                                                      % GU experiment #
OSU = [386];                % OSU experiment #
% OSU = [375 376 377 387 388 389 399 400 401 411 412 413 423 424 425 435 436 437];
% OSU = [378 379 380 390 391 392 402 403 404 414 415 416 426 427 428 438 439 440];
% OSU = [381 382 383 393 394 395 405 406 407 417 418 419 429 430 431 441 442 443];
other = 26;                                                                         % Saved reference case from other sources
gust_test = [3,6,11];                                                                   % Gust test
flap_test = 1:2;                                                                    % Flap test
tvfs_test = 14;                                                                   % TVFS test
hargen_test = {"tvfs_test-21"};                                                                % General harmonic test
% Select source as "NASA", "GU", "OSU", "other", "gust_test", "flap_test", "tvfs_test or "hargen_test"
INPUTS.source = "NASA";
% Select model ("BL","BLT","BLO","BLS","BLSLR","BLG" are available)
INPUTS.model = "BL";
% Select gust indicial response model as "K", "CFD", "BR-C" or "BR-F" (scaled Kussner, CFD from Leishman's book, or Berci and' Righi's)
INPUTS.gust_ind = "K";
% Select error standard as 'RMSE' or 'NRMSE'
INPUTS.error_standard = 'RMSE';
% Plotting separate figures option
INPUTS.plot_sep_figs = 0; 
INPUTS.save_folder = pwd;
INPUTS.figures_extension = '.pdf';

%% Call solver
call_solver;