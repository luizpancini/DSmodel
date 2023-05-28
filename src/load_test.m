function [params,data] = load_test(test)

% Standard values for b and a_inf if not input
if length(test) == 5
    b = 0.61/2;     % Same as NASA's experiment (see McAlister et al. (1982) - page 14)
    a_inf = 340;    % Assumed MSL conditions
elseif length(test) == 6
    b = test{6};
    a_inf = 340;
elseif length(test) == 7
    b = test{6};
    a_inf = test{7};
else
    error('Set test as a cell with at least 5 and at most 7 entries')
end

% Set remaining test conditions
M = test{1}; k = test{2}; a_0 = test{3}; a_1 = test{4}; airfoil = test{5}; 

% Flow properties
U = M*a_inf;        % Airspeed
beta = sqrt(1-M^2); % Prandtl-Glauert compressibility factor

% Set all flow and test condition variables into the params struct
params.M = M;
params.U = U;
params.b = b;
params.a_inf = a_inf;
params.beta = beta;
params.a_0 = a_0;
params.a_1 = a_1;
params.k = k;
params.airfoil = airfoil;

% Initialize data struct
data.authors = {[]};
        
end