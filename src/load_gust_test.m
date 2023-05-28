function [params,data] = load_gust_test(input)

% Load gust data for current input test case
gust_test_data = set_gust_test(input);

% Unpack current gust test conditions
airfoil = gust_test_data.airfoil;
ah = gust_test_data.ah;
a_0 = gust_test_data.a_0;
a_1 = gust_test_data.a_1;
a_inf = gust_test_data.a_inf;
b = gust_test_data.b;
beta = gust_test_data.beta;
k = gust_test_data.k;
M = gust_test_data.M;
U = gust_test_data.U;
lambda_g = gust_test_data.lambda_g;
gust_profile = gust_test_data.gust_profile;
t_init = gust_test_data.t_init;
tspan = gust_test_data.tspan;
if isfield(gust_test_data,'plot_increment')
    plot_increment = gust_test_data.plot_increment;
else
    plot_increment = true;
end
a_1h = 0; 
k_h = 0;

% Get gust profile variables
gust_options = gust_profile_vars(gust_profile,gust_test_data);

% Set parameters structure

params = variables2struct(struct(),airfoil,ah,a_0,a_1,a_inf,b,beta,k,M,U,lambda_g,gust_profile,gust_options,t_init,tspan,plot_increment,a_1h,k_h);

% Set experimental data structure
data = gust_test_data.data;

end