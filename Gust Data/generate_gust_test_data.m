clc
clear

% Set gust test conditions
gust_profile = "1-cos";
airfoil = "NACA0012";
M = 0.2;
nu = 1.845e-5;
a_inf = 295.07;
U = M*a_inf;
k = 0;
b = 1.22;
beta = sqrt(1-M^2);
a_0 = deg2rad(11);
a_1 = 0;
alpha_g = 4*pi/180;
wg_0 = U*atan(alpha_g);
H = 15/2;
semichords_travel_initialization = 50;
semichords_simulation = 40;
t_init = semichords_travel_initialization*b/U;
tf = semichords_simulation*b/U;
tspan = [0; tf];
f_sin = 1/50*U/b;

% Experimental data
authors = {"Mallik and Raveh (2020) - CFD", "Mallik and Raveh (2020) - ROM"};
filename1 = 'exp_data.mat'; filename2 = 'rom_data.mat';
load(filename1); load(filename2);
% data = gust_tests{1}.data;
data.authors = authors;
data.time_exp_cl = exp_data(:,1);
data.clt_exp = exp_data(:,2);
data.time_mod_cl = rom_data(:,1);
data.clt_mod = rom_data(:,2);

% Set all gust data on structure
gust_test = variables2struct(struct(),airfoil,a_0,a_1,a_inf,b,beta,k,M,U,gust_profile,t_init,tspan,H,wg_0,data);

% Prepare output as cell
gust_tests = {gust_test};

% Save data
filename_new = '1CG_H_15_A11.mat';
save(filename_new,'gust_tests');