clc
clear 

% Frames to be modified
data_list = 372:443;

% Loop
for data=data_list(:)'
    % Check filename
    filename = ['OSU_' num2str(data) '.mat'];
    if isfile(filename)
        % Load frame's data
        load(filename)
        % Set modifications
        a_0 = alpha_0;
        a_1 = delta_alpha;
        ah = -1/2;
        a_inf = 340;
        % Save
        save(filename,'authors','airfoil','b','a_inf','Re','M','U','beta','k','a_0','a_1','ah','time','alpha','cl','cd','cm','cn','cc','time_cycles','alpha_cycles','cl_cycles','cd_cycles','cm_cycles','cn_cycles','cc_cycles');
    end
end