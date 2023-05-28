clc
clear

% Set gust test conditions
gust_profile = "sharp-edge";
airfoil = "NACA0002";
M = 0.01;
a_inf = 340;
U = M*a_inf;
k = 0;
b = 0.5;
beta = sqrt(1-M^2);
a_0 = 0;
a_1 = 0;
alpha_g = 1*pi/180;
wg_0 = U*atan(alpha_g);
H = inf;
semichords_travel_initialization = 0;
semichords_simulation = 6;
t_init = semichords_travel_initialization*b/U;
tf = semichords_simulation*b/U;
tspan = [0; tf];
f_sin = nan;
lambda_g = -2;
plot_increment = false;

% Experimental data
authors = {"Leishman (1997)"};
filename1 = 'analytical_data.mat'; filename2 = 'mod_data.mat';
load(filename1); load(filename2);
% data = gust_tests{1}.data;
data.authors = authors;
data.time_exp_cl = analytical_data(:,1);
data.clt_exp = analytical_data(:,2);
data.time_mod_cl = nan;
data.clt_mod = nan;

% Set all gust data on structure
gust_test = variables2struct(struct(),airfoil,a_0,a_1,a_inf,b,beta,k,M,U,lambda_g,gust_profile,t_init,tspan,H,wg_0,data,plot_increment);

% Prepare output as cell
gust_tests = {gust_test};

% Save data
filename_new = 'SEG_inc_lambdan2.mat';
save(filename_new,'gust_tests');