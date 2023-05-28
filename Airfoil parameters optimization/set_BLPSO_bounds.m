function [LB,UB,LBmat,UBmat,Vmax,Vmaxmat] = set_BLPSO_bounds(N_p,varargin)

% Handle input
if ~isempty(varargin)
    airfoil = varargin{1};
end

switch airfoil
    case 'S809'
        B = [ pi/180*-1.5,  pi/180*-1.0     % alpha_0L
              pi/180*18.0,  pi/180*21.2     % alpha_ds0
              pi/180*14.7,  pi/180*15.3     % alpha_ss
               pi/180*7.5,  pi/180*10.0     % alpha1_0
                      0.4,          1.0     % gamma_LS
              pi/180*-4.0,   pi/180*1.0     % delta_alpha_0
               pi/180*0.0,   pi/180*4.0     % delta_alpha_1
                      0.0,          0.8     % delta_alpha_2
                      0.3,          1.5     % nu_1
                      0.5,          1.5     % nu_2
                    0.010,        0.014     % c_d0
                   -0.040,       -0.035     % c_m0
             180/pi*0.106, 180/pi*0.110     % c_n_alpha
               pi/180*0.0,   pi/180*1.5     % d_cc
               pi/180*0.0,   pi/180*2.0     % d_cm
                     0.00,         0.20     % E0
                     -0.1,          0.6     % E1
                      0.5,          2.0     % g_v
                    -0.02,        -0.01     % K0
                    -0.18,        -0.05     % K1
                    -0.04,         0.02     % K2
                     0.00,         0.06     % K3
                    0.001,        0.010     % r0
               pi/180*1.5,   pi/180*3.0     % S1
               pi/180*3.0,   pi/180*5.5     % S2
                      1.0,          3.5     % Ta
                      2.0,          6.0     % Tf0
                      3.5,          5.5     % TvL
                     0.20,         0.35     % Vm
                      2.0,          3.0     % Vn1
                      0.8,          1.6     % Vn2
                      0.8,          2.5     % z_cc
                      0.0,          2.0     % z_cm
                     0.25,         0.50     % gamma_TvL
                      2.0,          4.0     % kappa
                     0.00,         0.02     % df0_c
                     0.00,         0.02     % f0
                     0.50,         0.75     % fb
                      0.0,          1.0     % fSig1n
                      0.0,          1.0     % fSig1c
                      0.0,          2.0     % fSig2n
                      0.0,          1.0     % fS2n_ulpr
                      0.0,          2.5     % fS2n_dlpr
                    -0.25,         0.25     % fS1n_u
                    -0.50,         1.00     % fS1m_u
                    -0.25,         1.25     % fS1c_u
                    -0.35,          2.0     % fS1n_d
                    -0.75,         1.50     % fS1m_d
                    -0.35,          2.0     % fS1c_d
                      2.5,          5.0     % fS1n_ud
                      3.0,          5.0     % fS1m_ud
                      3.0,          6.0     % fS1c_ud
                    -0.25,          1.0     % fS2n_u
                      0.0,          4.0     % fS2m_u
                    -0.50,          1.0     % fS2c_u
                    -0.75,         0.50     % fS2n_d
                      0.5,          3.0     % fS2m_d
                    -0.50,          2.0     % fS2c_d
                       20,           80     % fSS1
                      1/4,            1     % fSS2
                      0.0,          5.0     % K1_f
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