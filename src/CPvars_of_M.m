function K0 = CPvars_of_M(M)

% Estimated from JOSE's CFD for the NACA0006
Mach_range = [   0.3;    0.5;   0.65;    0.7;  0.75;   0.8;   0.85;    0.9];
K0_range   = [-0.006; -0.005; -0.004; -0.003; 0.000; 0.005; -0.020; -0.190];

K0 = lininterp1(Mach_range,K0_range,M);

end