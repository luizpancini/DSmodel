function [alpha,c_n,c_m,c_c,c_l,c_d,t] = run_BLS(input,lin_reat,dt,n_cycles)

%% Load airfoil constants
[~,~,U,beta,M,k,a_0,a_1,b,~,~,c_n_alpha,S1,S2,alpha_ss,...
 alpha_ds0,Df,TvL,K0,K1,K2,~,kappa,Ta,Tf0,r0,E0,c_d0,alpha_0L,c_m0,f0,...
 fb,B1,B2,alpha1_0,alpha_min0,Tr,r_dcn_da_reat,r_dcm_da_reat,r_dcc_da_reat] = read_data_BLS(input);

%% Indicial parameters
[A1,A2,A3,b1,b2,b3,T_I,T_M] = read_indicial_params_BLS(M);

%% State space matrices
[A,B] = get_SS_matrices_BLS(U,b,M);

%% Initial conditions: {x0} = 0
N = 9; % Number of states
x0 = zeros(N,1);
y0 = [a_0; 0; 0; 0];
y0(5:20) = nan;

%% ODE solver
t_cycle = 2*pi*b/(U*k); % Time of one cycle
tf = t_cycle*n_cycles;  % Test time
tspan = [0 tf];
options = {'hlim',dt};
[t,~,y,~] = BLS_RKF45(tspan,x0,y0,A,B,U,M,b,a_0,a_1,k,beta,c_n_alpha,S1,S2,alpha1_0,alpha_ss,alpha_ds0,Df,TvL,K0,K1,K2,kappa,Ta,Tf0,r0,E0,c_d0,c_m0,f0,fb,B1,B2,alpha_0L,T_I,T_M,b1,b2,b3,A1,A2,A3,alpha_min0,Tr,r_dcn_da_reat,r_dcm_da_reat,r_dcc_da_reat,lin_reat,options);

%% Output variables
alpha = y(1,:);
c_n = y(2,:);
c_m = y(3,:);
c_c = y(4,:);
c_l = c_n.*cos(alpha)+c_c.*sin(alpha);
c_d = c_n.*sin(alpha)-c_c.*cos(alpha);

end