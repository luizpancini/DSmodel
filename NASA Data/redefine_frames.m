clc
clear 

% Frames to be modified
frame_list = 7019:31310;

% Loop
for frame=frame_list(:)'
    % Check filename
    filename = ['frame_' num2str(frame) '.mat'];
    if isfile(filename)
        % Load frame's data
        load(filename)
        % Set modifications
        a_0 = alpha_0;
        a_1 = delta_alpha;
        ah = -1/2;
        % Save
        save(filename,'authors','airfoil','a_inf','b','ah','U','M','beta','a_0','a_1','k','alpha_exp_cl','cl_exp','alpha_exp_cm','cm_exp','alpha_exp_cd','cd_exp','time_exp_cl','clt_exp','time_exp_cm','cmt_exp','time_exp_cd','cdt_exp');
    end
end