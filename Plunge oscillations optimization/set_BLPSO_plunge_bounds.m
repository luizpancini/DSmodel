function [D,LB,UB,LBmat,UBmat,Vmax,Vmaxmat] = set_BLPSO_plunge_bounds(N_p,varargin)

% Handle input
if ~isempty(varargin)
    airfoil = varargin{1};
end

switch airfoil
    case 'Vertol-23010'
        B = [  pi/180*0.5,   pi/180*0.7     % alpha_0L
              pi/180*16.0,  pi/180*21.0     % alpha_ds0
              pi/180*13.2,  pi/180*13.6     % alpha_ss
              pi/180*12.0,  pi/180*15.0     % alpha1_0
                      0.1,          1.0     % gamma_LS
               pi/180*0.5,   pi/180*4.0     % delta_alpha_0
               pi/180*0.5,   pi/180*5.0     % delta_alpha_1
                      0.0,          0.0     % delta_alpha_2
                      0.3,          1.5     % nu_1
                      0.5,          1.5     % nu_2
                    0.007,        0.012     % c_d0
                   -0.015,       -0.005     % c_m0
                2*pi/0.93,    2*pi/0.88     % c_n_alpha
               pi/180*0.0,   pi/180*0.0     % d_cc
               pi/180*0.0,   pi/180*0.0     % d_cm
                     0.10,         0.40     % E0
                     -0.1,          0.5     % E1
                      0.5,          2.0     % g_v
                    -0.03,         0.02     % K0
                    -0.18,        -0.05     % K1
                    -0.03,         0.08     % K2
                     0.00,         0.30     % K3
                    0.005,        0.030     % r0
               pi/180*1.0,   pi/180*3.0     % S1
               pi/180*1.0,   pi/180*4.0     % S2
                      0.5,          2.5     % Ta
                      1.5,          5.0     % Tf0
                      3.0,          8.0     % TvL
                     0.20,         0.45     % Vm
                      0.5,          1.2     % Vn1
                      0.1,          0.5     % Vn2
                      0.9,          1.1     % z_cc
                      0.0,          1.5     % z_cm
                     0.25,         0.50     % gamma_TvL
                      1.5,          4.0     % kappa
                     0.00,         0.02     % df0_c
                     0.00,         0.02     % f0
                     0.65,         0.85     % fb
                      0.0,          2.0     % fSig1n
                      0.0,          0.0     % fSig1c
                      0.0,          2.0     % fSig2n
                      0.0,          1.0     % fS2n_ulpr
                      0.0,          2.5     % fS2n_dlpr
                    -0.25,         0.25     % fS1n_u
                    -0.50,         1.00     % fS1m_u
                    -0.75,        -0.75     % fS1c_u
                    -0.35,          2.0     % fS1n_d
                    -0.75,         1.50     % fS1m_d
                      1.0,          1.0     % fS1c_d
                      1.0,          5.0     % fS1n_ud
                      1.0,          6.0     % fS1m_ud
                      2.0,          2.0     % fS1c_ud
                    -0.35,          1.0     % fS2n_u
                      0.0,          2.0     % fS2m_u
                      2.0,          2.0     % fS2c_u
                    -0.75,         0.50     % fS2n_d
                      0.3,          2.0     % fS2m_d
                      2.0,          2.0     % fS2c_d
                       10,           80     % fSS1
                      1/4,            1     % fSS2
                      0.0,          2.0     % K1_f
                                       ];

end

% Number of optimization variables
D = size(B,1);

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