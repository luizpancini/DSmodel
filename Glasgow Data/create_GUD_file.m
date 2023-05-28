clear

%% Select files
GUD_range = [11000011:10:11000861,11004131:10:11004471,11011962:10:11012862,11012891:10:11014461]; % NACA0012
GUD_range = [12000411:10:12000831,12001111:10:12001721,12010702:10:12010772]; % NACA0015-s
GUD_range = [06000011:10:06000081,06011001:10:06014261,06014272:10:06014462,06014391:10:06014421,06019001:10:06019561]; % NACA0018
GUD_range = [05000011:10:05000081,05011001:10:05019561]; % NACA0015
GUD_range = [02000101,02010021:10:02010931]; % NACA23012A

for GUD = GUD_range(:)'
    filename = sprintf('%08d.dat',[GUD]); 
    filename_coeffs = sprintf('%08d_coeffs.dat',[GUD]);   
    if isfile(filename_coeffs) && isfile(filename)
        % Get test conditions data
        fid = fopen(filename);
        RIB = fscanf(fid,'%g',32); % read in first 32 elements of RIB
        fclose('all');
        M = RIB(21);
        k = RIB(22);
        U = RIB(23);
        ah = -1/2;
        a_inf = U/M;
        beta = sqrt(1-M^2);
        authors = {'Green & Giuni (2017)'};
        airfoil_code = filename(1:2);
        switch airfoil_code
            case '02'
                airfoil = 'NACA23012A';
                b = 0.55/2;
            case '05'
                airfoil = 'NACA0015';
                b = 0.55/2;
            case '06'
                airfoil = 'NACA0018';
                b = 0.55/2;
            case '11'
                airfoil = 'NACA0012-GU';
                b = 0.55/2;
            case '12'
                airfoil = 'NACA0015-s';
                b = 0.55/4;
        end
        % Nominal alpha_0 and delta_alpha are RIB(8) and RIB(9), but we
        % shall use the ones gotten from the coefficients data
        % Get coefficients data
        GUD_data = importdata(filename_coeffs); 
        alpha = GUD_data.data(:,2);
        cn = GUD_data.data(:,3);
        cc = GUD_data.data(:,4);
        cm = GUD_data.data(:,5);
        cl = cn.*cosd(alpha)+cc.*sind(alpha);
        cd = cn.*sind(alpha)-cc.*cosd(alpha);
        a_0 = (max(alpha)+min(alpha))/2*pi/180;
        a_1 = (max(alpha)-min(alpha))/2*pi/180;
        if k > 0
            time = GUD_data.data(:,1);      % cycle angle (0 -- 2*pi)
            time = time/(2*pi)*360-90;      % transform to deg (-90 -- 270)
        else
            k = 1e-3;
            [~,ind] = max(alpha); 
            time = real(asind((alpha-a_0*180/pi)/(a_1*180/pi)));
            time(ind+1:end) = 180-time(ind+1:end); 
            if time(65) == 90, time(65) = 89.99; end
        end
        % Save data
        filename_save = sprintf('GUD_%08d.mat', GUD); 
        save(filename_save,'authors','airfoil','b','U','beta','a_inf','M','k','a_0','a_1','ah','time','alpha','cn','cc','cm','cl','cd');
%         pause(0.5)
    end
end