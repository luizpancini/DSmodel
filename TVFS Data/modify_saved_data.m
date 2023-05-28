clc
clear

% Set cases to modify
cases = 1:23;

for case_now=cases(:)'
    % Load file
    filename = "tvfs_case_" + num2str(case_now,'%02d') + ".mat";
    load(filename);
    % Set modifications
    ah = -1/2;
    % Save file
    save(filename,'authors','M','U','b','a_inf','beta','a_0','a_1','k','U_0','U_1','k_U','psi_a','airfoil','lambda_U','ah','time_exp_cl','clt_exp','time_mod_cl','clt_mod','time_exp_cm','cmt_exp','time_mod_cm','cmt_mod');
end
