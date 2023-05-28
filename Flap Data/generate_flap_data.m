clc
clear

% Experimental setup
flap_case = "01";
authors = {"Tijdeman & Schippers (1973) - exp.","Leishman (2006) - mod."};
airfoil = "NACA64A006";
M = 0.5;
d_0 = 0;
d_1 = 2.5*pi/180;
kf = 0.098;
a_inf = 340;
U = M*a_inf;
beta = sqrt(1-M^2);
dh = 1/2;
ah = -1/2;  % This is unknown, but likely
b = 0.3;    % This is unknown, but does not influence the incompressible solution

% Available coefficients data
coefs = ["cn"; "cm"; "ch"];

% Load data points
for i=1:length(coefs)
    % Load
    filename_exp = "exp_" + coefs(i) + "_" + flap_case + ".mat";
    filename_mod = "mod_" + coefs(i) + "_" + flap_case + ".mat";
    load(filename_exp); 
    load(filename_mod);
end

% Get constants to multiply data (scaled plots)
f_n = pi;
f_m = -pi/2;
f_h = -pi/2;

% Set data arrays
delta_exp_cn = exp_cn(:,1);
delta_mod_cn = mod_cn(:,1);
cn_exp = f_n*exp_cn(:,2);
cn_mod = f_n*mod_cn(:,2);
delta_exp_cm = exp_cm(:,1);
delta_mod_cm = mod_cm(:,1);
cm_exp = f_m*exp_cm(:,2);
cm_mod = f_m*mod_cm(:,2);
delta_exp_ch = exp_ch(:,1);
delta_mod_ch = mod_ch(:,1);
ch_exp = f_h*exp_ch(:,2);
ch_mod = f_h*mod_ch(:,2);

% Save file
new_filename = "flap_case_" + flap_case + ".mat";
save(new_filename,'authors','airfoil','b','a_inf','M','U','beta','ah','dh','kf','d_0','d_1','delta_exp_cn','delta_exp_cm','delta_exp_ch','cn_exp','cm_exp','ch_exp','delta_mod_cn','delta_mod_cm','delta_mod_ch','cn_mod','cm_mod','ch_mod');
