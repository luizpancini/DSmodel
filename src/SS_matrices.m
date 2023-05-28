function params = SS_matrices(params,model,gust_ind)

%% State-space matrices for each model
% Unpack
U = params.U;
b = params.b;
beta = params.beta;
switch model
    case {"BL","BLgust","BLtvfs"}
        % Unpack
        b1 = params.b1;
        b2 = params.b2;
        b3 = params.b3;
        b4 = params.b4;
        b5 = params.b5;
        K_a = params.K_a;
        K_q = params.K_q;
        K_M = params.K_M;
        K_aM = params.K_aM;
        K_qM = params.K_qM;
        K_MM = params.K_MM;
        T_I = params.T_I;
        % System matrix (diagonal only)
        params.A = [-U/b*beta^2*b1; -U/b*beta^2*b2; -1/(K_a*T_I); -1/(K_q*T_I); -1/(b3*K_aM*T_I); -1/(b4*K_aM*T_I); -b5*U/b*beta^2; -1/(K_qM*T_I); -1/(K_M*T_I); -1/(b3*K_MM*T_I); -1/(b4*K_MM*T_I);];
        % Input matrix
        params.B = [1 1/2; 1 1/2; 1 0; 0 1; 1 0; 1 0; 0 1; 0 1];
    case {"BLO","BLG","BLOgust"}
        % Unpack
        b1 = params.b1;
        b2 = params.b2;
        b3 = params.b3;
        b4 = params.b4;
        b5 = params.b5;
        K_a = params.K_a;
        K_q = params.K_q;
        K_aM = params.K_aM;
        K_qM = params.K_qM;
        T_I = params.T_I;
        % System matrix (diagonal only)
        params.A = [-U/b*beta^2*b1; -U/b*beta^2*b2; -1/(K_a*T_I); -1/(K_q*T_I); -1/(b3*K_aM*T_I); -1/(b4*K_aM*T_I); -b5*U/b*beta^2; -1/(K_qM*T_I)];
        % Input matrix
        params.B = [1 1/2; 1 1/2; 1 0; 0 1; 1 0; 1 0; 0 1; 0 1];    
    case {"BLT","BLTgust","BLTflap","BLTtvfs"}
        % Unpack
        b1 = params.b1;
        b2 = params.b2;
        c2 = params.c2;
        c4 = params.c4;
        % System matrix (diagonal only)
        params.A = -U/b*beta^2*[b1; b2];
        % Input matrix
        params.B = [1 1/2; 1 1/2];
%         params.B = -U/b*[c2*c4*U/b, c2+c4]; % For thin-airfoil theory potential flow states
    case {"BLS","BLSLR"}
        % Unpack
        b1 = params.b1;
        b2 = params.b2;
        b3 = params.b3;
        T_I = params.T_I;
        T_M = params.T_M;
        % System matrix (diagonal only)
        params.A = -U/b*[b1; b2; b3; 1/T_I; 1*T_M; 1/T_I];
        % Input matrix
        params.B = [1 1/2; 1 1/2; 1 1/2; 1 1/4; 0 1/2; 1/2 7/24];
end

%% Gust modeling state-space matrices (this system is equivalent to Leishman's eqs. 8.195-196)
if isfield(params,'AG')
    params.gust_ind = gust_ind;
    % Circulatory-only matrices
    AG = params.AG; bG = params.bG;     % Exponential approximation coefficients
    params.Ag = -U/b*beta^2*diag(bG);   % System matrix
    params.Bg = ones(length(AG),1);     % Input matrix
    params.Cg = U/b*beta^2*(AG.*bG)';   % Output matrix (circulatory normal coefficient)
    params.Cg2 = [];                    % Output matrix (inertial normal coefficient)
    % Inertial contributions from Berci and Righi's model
    if gust_ind == "BR-F"
        params.Ag = -U/b*beta^2*diag([bG(1:7); bG(7); bG(8:end)]); params.Ag(7,8) = -U/b*beta^2*params.omega1I; params.Ag(8,7) = U/b*beta^2*params.omega1I;
        params.Bg = ones(11,1); params.Bg(8) = 0;
        params.Cg = U/b*beta^2*(AG(1:6).*bG(1:6))';
        params.Cg2 = -U/b*beta^2*[AG(7)*bG(7), AG(7)*params.omega1I, AG(8)*bG(8), AG(9)*bG(9), AG(10)*bG(10)];
    end
    % Modifications for general harmonic model with semi-empirical gust-convection c_n
    if contains(model,'gen') && isfield(params,'g1')
        params.bG = [params.g1; params.g2; params.bG];
    end
end

end