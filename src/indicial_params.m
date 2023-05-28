function params = indicial_params(params,model,gust_ind)

% Unpack
M = params.M;
b = params.b;
a_inf = params.a_inf;
beta = params.beta;
airfoil = params.airfoil;

%% Indicial parameters according to each model
switch model
    case {"BL","BLgust","BLtvfs"}
        %% Circulatory indicial parameters
        % Nominal values
        params.A1 = 0.3;
        params.A2 = 0.7;
        params.A3 = 1.5;
        params.A4 = -0.5;
        params.b1 = 0.14;
        params.b2 = 0.53;
        params.b3 = 0.25;
        params.b4 = 0.1;
        params.b5 = 5.0;
        % Adjustments for each airfoil
        switch airfoil
            case {'NACA0012','NACA0012-GU','NACA0015','NACA0015-s','NACA0018','NACA23012A'}
                params.b1 = params.b1*1.5;
                params.b2 = params.b2*1.0;
            case 'AMES-01'
                params.b1 = params.b1*3.0;
                params.b2 = params.b2*0.6;
            case 'NLR-7301'
                params.b1 = params.b1*2.5;
                params.b2 = params.b2*0.8;
            case 'S809'
                params.b1 = params.b1*1.0;
                params.b2 = params.b2*1.0;
            case 'Vertol-23010'
                params.b1 = params.b1*1.0;
                params.b2 = params.b2*1.0;    
        end
        % Modifications for time-varying freestream tests - Use Jose's approximation for the circulatory indicial function
        if contains(model,"tvfs")
            params.A1 = 0.3493;
            params.A2 = 0.6507;
            params.b1 = 0.0984;
            params.b2 = 0.7759;
        % Modifications for flap tests - Use Wagner's approximation for the circulatory indicial function    
        elseif contains(model,"flap")
            params.A1 = 0.165;
            params.A2 = 0.335;
            params.b1 = 0.0455;
            params.b2 = 0.3;    
        end
        %% Impulsive indicial parameters
        % Nominal values
        params.K_a = 1/(1-M+pi*beta*M^2*(params.A1*params.b1+params.A2*params.b2));
        params.K_q = 1/(1-M+2*pi*beta*M^2*(params.A1*params.b1+params.A2*params.b2));
        params.K_M = 1/(1-M+pi/beta*M^2*(params.A1*params.b1+params.A2*params.b2));
        params.K_aM = (params.A3*params.b4+params.A4*params.b3)/(params.b3*params.b4*(1-M));
        params.K_qM = 7/(15*(1-M)+3*pi*beta*M^2*params.b5);
        params.K_MM = params.K_aM;
        params.T_I = 2*b/a_inf;
        % Leishman and Nguyen (1989) say they have reduced the K constants by 25% to
        % match experimental data (see Conclusions of their paper). Here we apply
        % that only for M > 0.07, and only for c_n related impulsive parameters.
        if M > 0.07
            fac = 0.75;
        else
            fac = 1-0.25*(M/0.07)^2;
        end
        params.K_a = params.K_a*fac;
        params.K_q = params.K_q*fac;
        params.K_M = params.K_M*fac;
        % Adjustments for each airfoil
        switch airfoil
            case 'NACA0012'
                params.K_aM = params.K_aM*1.0;
                params.K_qM = params.K_qM*1.0;
                params.K_MM = params.K_MM*1.0;
            case 'AMES-01'
                params.K_aM = params.K_aM*1.25;
                params.K_qM = params.K_qM*1.25;
                params.K_MM = params.K_MM*1.25;
            case 'NLR-7301'
                params.K_aM = params.K_aM*1.0;
                params.K_qM = params.K_qM*1.0;
                params.K_MM = params.K_MM*1.0;
            case 'S809'
                params.K_aM = params.K_aM*2.0;
                params.K_qM = params.K_qM*2.0;
                params.K_MM = params.K_MM*2.0;
            case 'Vertol-23010'
                params.K_aM = params.K_aM*1.15;
                params.K_qM = params.K_qM*1.15;
                params.K_MM = params.K_MM*1.15;    
        end
        % Convecting gust parameters
        if contains(model,["gust"])
            params = convecting_gust_params(params,gust_ind);
        end
    case {"BLT","BLTgust","BLTflap","BLTtvfs"}
        % Wagner function constants
        params.c1 = 0.165;
        params.c2 = 0.0455;
        params.c3 = 0.335;
        params.c4 = 0.3;
        % Circulatory indicial parameters
        % Nominal values
        params.A1 = 0.3;
        params.A2 = 0.7;
        params.b1 = 0.14;
        params.b2 = 0.53;
        % Adjustments for each airfoil
        switch airfoil
            case 'NACA0012'
                params.b1 = params.b1*2.0;
                params.b2 = params.b2*1.0;
            case 'AMES-01'
                params.b1 = params.b1*3.0;
                params.b2 = params.b2*0.6;
            case 'NLR-7301'
                params.b1 = params.b1*2.5;
                params.b2 = params.b2*0.8;
            case 'S809'
                params.b1 = params.b1*1.0;
                params.b2 = params.b2*1.0;
            case 'NACA64A006'
                params.b1 = params.b1*2.5;
                params.b2 = params.b2*1.25;
        end
        % Modifications for time-varying freestream and flap tests - Use Wagner's approximation for the circulatory indicial function 
        if contains(model,"tvfs")
            params.A1 = 0.165;
            params.A2 = 0.335;
            params.b1 = 0.0455;
            params.b2 = 0.3;
        elseif contains(model,["flap"])
            params.A1 = 0.165;
            params.A2 = 0.335;
            params.b1 = 0.0455;
            params.b2 = 0.3;
        end
    case {"BLS","BLSLR"}
        params.A1 = 0.165;
        params.A2 = 0.335;
        params.A3 = 0.5;
        params.b1 = 0.05;
        params.b2 = 0.222;
        params.b3 = 0.8/M;
        params.T_M = M/2;
        params.T_I = M*(1+3*M)/2;
    case {"BLO","BLG","BLOgust"}
        % Circulatory indicial parameters
        params.A1 = 0.3;
        params.A2 = 0.7;
        params.A3 = 1.5;
        params.A4 = -0.5;
        params.b1 = 0.14;
        params.b2 = 0.53;
        params.b3 = 0.25;
        params.b4 = 0.1;
        params.b5 = 0.5;
        % Impulsive indicial parameters
        params.K_a = 0.75/(1-M+pi*beta*M^2*(params.A1*params.b1+params.A2*params.b2));
        params.K_q = 0.75/(1-M+2*pi*beta*M^2*(params.A1*params.b1+params.A2*params.b2));
        params.K_aM = (params.A3*params.b4+params.A4*params.b3)/(params.b3*params.b4*(1-M));
        params.K_qM = 7/(15*(1-M)+3*pi*beta*M^2*params.b5);
        params.T_I = 2*b/a_inf;
    case "BLhargen"
        %% Circulatory indicial parameters
        % Nominal values
        params.A1 = 0.3;
        params.A2 = 0.7;
        params.A3 = 1.5;
        params.A4 = -0.5;
        params.b1 = 0.14;
        params.b2 = 0.53;
        params.b3 = 0.25;
        params.b4 = 0.1;
        params.b5 = 5.0;
        % Adjustments for each airfoil
        switch airfoil
            case 'NACA0012'
                params.b1 = params.b1*1.5;
                params.b2 = params.b2*1.0;
            case 'AMES-01'
                params.b1 = params.b1*3.0;
                params.b2 = params.b2*0.6;
            case 'NLR-7301'
                params.b1 = params.b1*2.5;
                params.b2 = params.b2*0.8;
            case 'S809'
                params.b1 = params.b1*1.0;
                params.b2 = params.b2*1.0;
        end
        % Modifications for time-varying freestream tests - Use Jose's approximation for the circulatory indicial function
        if params.mod_ind_params
            params.A1 = 0.3493;
            params.A2 = 0.6507;
            params.b1 = 0.0984;
            params.b2 = 0.7759;
        end
        % Flap motion parameters - Use Wagner's approximation for the circulatory indicial function
        params.A1W = 0.165;
        params.A2W = 0.335;
        params.b1W = 0.0455;
        params.b2W = 0.3;
        % Convecting gust parameters
        params = convecting_gust_params(params,gust_ind);
    case "BLThargen"
        %% Circulatory indicial parameters
        % Nominal values
        params.A1 = 0.3;
        params.A2 = 0.7;
        params.b1 = 0.14;
        params.b2 = 0.53;
        % Adjustments for each airfoil
        switch airfoil
            case 'NACA0012'
                params.b1 = params.b1*2.5;
                params.b2 = params.b2*1.0;
            case 'AMES-01'
                params.b1 = params.b1*3.0;
                params.b2 = params.b2*0.6;
            case 'NLR-7301'
                params.b1 = params.b1*2.5;
                params.b2 = params.b2*0.8;
            case 'S809'
                params.b1 = params.b1*1.0;
                params.b2 = params.b2*1.0;
            case 'NACA64A006'
                params.b1 = params.b1*2.5;
                params.b2 = params.b2*1.25;
        end
        % Modifications for time-varying freestream tests - Use Jose's approximation for the circulatory indicial function
        if params.mod_ind_params
            params.A1 = 0.3493;
            params.A2 = 0.6507;
            params.b1 = 0.0984;
            params.b2 = 0.7759;
        end
        % Flap motion parameters - Use Wagner's approximation for the circulatory indicial function
        params.A1W = 0.165;
        params.A2W = 0.335;
        params.b1W = 0.0455;
        params.b2W = 0.3;
end

%% Gust model indicial parameters
if contains(model,'BLT')
    switch gust_ind
        case "K"    % Jones' approximation to Wagner's function - this resumes to Kussner's function in the case of the stationary gust (lambda_g=1)
            params.N_g = 2;
            params.AG = [0.165; 0.335];
            params.bG = [0.0455; 0.3];
        otherwise
            error("Select 'K' as the gust indicial model when using the 'BLT' model")
    end
else
    switch gust_ind
        case "K"    % Jones' approximation to Kussner's function
            params.N_g = 2;
            params.AG = [0.5; 0.5];
            params.bG = [0.13; 1.0];
        case "CFD"  % Leishman's book CFD (Table 8.3)
            params.N_g = 2;
            params.AG = [0.67; 0.33];
            params.bG = [0.1753; 1.637];
        case "BR-C"   % Berci and Righi's with circulatory part only (An enhanced analytical method for the subsonic indicial lift of two-dimensional aerofoils – with numerical cross-validation)
            params.N_g = 6;
            params.AG = [0.3694; 0.0550; 0.2654; 0.1829; 0.0861; 0.0412];
            params.bG = [0.3733; 0.0179; 0.1096; 1.6003; 10.428; 170.93];
        case "BR-F"   % Berci and Righi's with full loads: circulatory and inertial (An enhanced analytical method for the subsonic indicial lift of two-dimensional aerofoils – with numerical cross-validation)
            params.N_g = 10;
            % Circulatory
            A1 = 0.3694; A2 = 0.0550; A3 = 0.2654; A4 = 0.1829; A5 = 0.0861; A6 = 0.0412;
            b1 = 0.3733; b2 = 0.0179; b3 = 0.1096; b4 = 1.6003; b5 = 10.428; b6 = 170.93;
            % Inertial
            c_n_alpha = params.c_n_alpha;
            A7 = -c_n_alpha*(A5+A6+A4); A8 = c_n_alpha*A6; A9 = c_n_alpha*A5; A10 = c_n_alpha*A4;
            b7 = 1/A7*(c_n_alpha*(A1*b1+A2*b2+A3*b3)-2/(beta^2*sqrt(M))); b8 = b6; b9 = b5; b10 = b4;
            omega1I = sqrt(b7^2-(c_n_alpha/A7*(A1*b1^2+A2*b2^2+A3*b3^2)));
            params.AG = [A1; A2; A3; A4; A5; A6; A7; A8; A9; A10];
            params.bG = [b1; b2; b3; b4; b5; b6; b7; b8; b9; b10];
            params.omega1I = omega1I;
    end
end

%% Flap model indicial parameters (Theodorsen's flap coefficients)
if contains(model,["flap","gen"])
    ah = params.ah;     % Semichord-normalized pitch axis position after midchord
    dh = params.dh;     % Semichord-normalized hinge position after midchord
    params.T1 = -1/3*sqrt(1-dh^2)*(2+dh^2)+dh*acos(dh);
    params.T2 = dh*(1-dh^2)-sqrt(1-dh^2)*(1+dh^2)*acos(dh)+dh*acos(dh)^2;
    params.T3 = -(1/8+dh^2)*acos(dh)^2+1/4*dh*sqrt(1-dh^2)*acos(dh)*(7+2*dh^2)-1/8*(1-dh^2)*(5*dh^2+4);
    params.T4 = -acos(dh)+dh*sqrt(1-dh^2);
    params.T5 = -(1-dh^2)-acos(dh)^2+2*dh*sqrt(1-dh^2)*acos(dh);
    params.T6 = params.T2;
    params.T7 = -(1/8+dh^2)*acos(dh)+1/8*dh*sqrt(1-dh^2)*(7+2*dh^2);
    params.T8 = -1/3*sqrt(1-dh^2)*(2*dh^2+1)+dh*acos(dh);
    params.T9 = 1/2*(1/3*sqrt(1-dh^2)^3+ah*params.T4);
    params.T10 = sqrt(1-dh^2)+acos(dh);
    params.T11 = acos(dh)*(1-2*dh)+sqrt(1-dh^2)*(2-dh);
    params.T12 = sqrt(1-dh^2)*(2+dh)-acos(dh)*(2*dh+1);
    params.T13 = 1/2*(-params.T7-(dh-ah)*params.T1);
end

end