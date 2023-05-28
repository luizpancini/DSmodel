function [LB,UB,LBmat,UBmat,Vmax,Vmaxmat] = set_BLPSO_cgust_bounds(N_p,varargin)

% Handle input
if ~isempty(varargin)
    airfoil = varargin{1};
end

switch airfoil
    case 'NACA0006'
        B = [  0,   30     % G @ lambda = 0.5
               0,   30     % G @ lambda = 0.6
               0,   30     % G @ lambda = 0.7
               0,   30     % G @ lambda = 0.75
             -20,    0     % G @ lambda = 1.5
             -30,    0     % G @ lambda = 2.0
             0.1,  5.0     % g1 @ lambda = 0.5
             0.1,  5.0     % g1 @ lambda = 0.6
             0.1,  5.0     % g1 @ lambda = 0.7
             0.1,  5.0     % g1 @ lambda = 0.75
             0.1,  5.0     % g1 @ lambda = 1.5
             0.1,  5.0     % g1 @ lambda = 2.0
             0.1,  5.0     % g2 @ lambda = 0.5
             0.1,  5.0     % g2 @ lambda = 0.6
             0.1,  5.0     % g2 @ lambda = 0.7
             0.1,  5.0     % g2 @ lambda = 0.75
             0.1,  5.0     % g2 @ lambda = 1.5
             0.1,  5.0     % g2 @ lambda = 2.0
                      ];

end

% Lower and upper bounds
LB = B(:,1)';
UB = B(:,2)';

% Check inconsistent bounds
if any(LB > UB)
    error('There are inconsistent bounds')
end

% Lower and upper bounds repeated over rows
LBmat = repmat(LB,N_p,1);
UBmat = repmat(UB,N_p,1);

% Maximum velocities (default is a quarter of the domain of the dimension)
Vmax = (UB-LB)/4;
Vmaxmat = repmat(Vmax,N_p,1);

end