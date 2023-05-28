clc
clear

% Load file
filename = "1CG_H_15_A11.mat";
load(filename,'gust_tests');

% Set modifications
load('analytical_data.mat');
gust_tests{1}.ah = -1/2;

% Save data
save(filename,'gust_tests');
