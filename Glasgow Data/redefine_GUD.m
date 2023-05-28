clc
clear 

% Frames to be modified
GUD_list = 11012891:10:11014461;

% Loop
for GUD=GUD_list(:)'
    % Check filename
    filename = ['GUD_' num2str(GUD) '.mat'];
    if isfile(filename)
        % Load frame's data
        load(filename)
        % Set modifications
        airfoil = 'NACA0012-GU';
        % Save
        save(filename,'authors','airfoil','b','U','beta','a_inf','M','k','a_0','a_1','ah','time','alpha','cn','cc','cm','cl','cd');
    end
end
