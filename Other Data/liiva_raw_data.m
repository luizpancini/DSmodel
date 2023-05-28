clc
clear
close all

%% Case selection
case_now = 3;
switch case_now
    case 1
        U = 426.1*.3048; k = 0.068; a_0 = 12.25; Delta_h = 0.306;
        cn_amp = [1.376; 0.116; 0.060; 0.021; 0.013; 0.011; 0.005; 0.005; 0.005; 0.005];
        cn_psi = [    0;   306;   170;   318;    64;   237;    58;   159;   300;    74];
    case 2
        U = 428.1*.3048; k = 0.121; a_0 = 12.46; Delta_h = 0.306;
        cn_amp = [1.279; 0.275; 0.074; 0.028; 0.020; 0.011; 0.003; 0.003; 0.003; 0.003];
        cn_psi = [    0;    305;   94;   156;    282;   142;    20;   315;   36;   199];
    case 3 % Page 300 (TP 4041.1)
        U = 435.5*.3048; k = 0.116; a_0 = 12.45; Delta_h = 0.472;
        cn_amp = [1.179; 0.341; 0.125; 0.039; 0.013; 0.019; 0.002; 0.011; 0.008; 0.008];
        cn_psi = [    0;   314;   130;   248;   348;   162;   331;   307;   137;   322];
        % Digitized reference data from volume 1
        load('CN_Vertol23010_M_04_k_0116_a0_12_a1h_3_exp.mat');
        alpha_exp_cn = CN_Vertol23010_M_04_k_0116_a0_12_a1h_3_exp(:,1);
        cn_exp = CN_Vertol23010_M_04_k_0116_a0_12_a1h_3_exp(:,2);
        % d = 25 to 35 appears to be close to the correct value for the phase in this case
    case 4 % Page 203 (TP 1102.4) - Pitching oscillation at k=0.355
        U = 444.3*.3048; k = 0.355; a_0 = 12.45; a_1 = 5.65;
        cn_amp = [1.081; 0.539; 0.062; 0.009; 0.015; 0.023; 0.009; 0.006; 0.006; 0.005];
        cn_psi = [    0;     9;   194;    49;   176;     2;   224;    90;   318;   41];
        % Digitized reference data from volume 1
        load('CN_Vertol23010_M_04_k_0355_a0_123_a1_5_exp.mat');
        alpha_exp_cn = CN_Vertol23010_M_04_k_0355_a0_123_a1_5_exp(:,1);
        cn_exp = CN_Vertol23010_M_04_k_0355_a0_123_a1_5_exp(:,2);
        % d = 170 appears to be close to the correct value for the phase in this case
end

%%
c = 6.38*0.0254; b = c/2;
om = k*U/b; % rad/s
T = 2*pi/om; % period
if ~exist('a_1','var')
    a_1 = atand(Delta_h*k);
end
t = T*linspace(0,1,1e4)';
harmonics = 1:9;
cn = zeros(length(t),1);
for i=1:9
    cn = cn + cn_amp(i+1)*sin(harmonics(i)*om*t+cn_psi(i+1));
end
cn = cn+cn_amp(1);
d = 30;
t_deg = t/T*360+d;
alpha_eq = a_0+a_1*sind(t_deg);
% figure;plot(t_deg,cn,'k-');grid; %xlim([-90,270]); ax = gca; ax.XTick = -90:90:270;
figure;plot(alpha_eq,cn,'k-'); hold on; grid; 
if exist('alpha_exp_cn','var')
    plot(alpha_exp_cn,cn_exp,'ko'); grid
end