clc
clear

% Data to be modified
data_list = 1:2;

% Loop
for data=data_list(:)'
    % Check filename
    filename = "flap_case_" + num2str(data,'%02d') + ".mat";
    if isfile(filename)
        % Load frame's data
        load(filename)
        % Set modifications
        a_0 = 0;
        a_1 = 0;
        k = 0;
        % Save
        save(filename,'authors','airfoil','b','a_inf','M','U','beta','ah','dh','k_f','d_0','d_1','a_0','a_1','k','delta_exp_cn','delta_exp_cm','delta_exp_ch','cn_exp','cm_exp','ch_exp','delta_mod_cn','delta_mod_cm','delta_mod_ch','cn_mod','cm_mod','ch_mod');
    end
end
