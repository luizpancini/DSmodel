clc
clear
% close all

% Constants
alpha1_0 = 8.2*pi/180;
alpha2_0 = 17.0*pi/180;
f01 = 1; f02 = 0; f03 = 0.02;
fb1 = -0.04995; fb2 = 2.8844; fb3 = 4.08e4;
S1 = 12.066;
S2 = log((f01+fb1*exp(S1*alpha1_0)-f02)/fb2)/alpha1_0;
S3 = log((f02+fb2*exp(S2*alpha2_0)-f03)/fb3)/alpha2_0;
K0 = -3e-3;
K1 = -1e-3;
K2 = -2.5e-2;
kappa = 6;

alpha = pi/180*(0:0.1:40);
f = nan(length(alpha),1);
dCP = f;
f_a2 = f02 + fb2*exp(S2*abs(alpha2_0));
for i=1:length(alpha)
    % f
    if abs(alpha(i)) <= alpha1_0
        f(i) = f01 + fb1*exp(S1*abs(alpha(i)));
    elseif abs(alpha) >= alpha2_0
        f(i) = f03 + fb3*exp(S3*abs(alpha(i)));
    else
        f(i) = f02 + fb2*exp(S2*abs(alpha(i)));
    end   
%     if f(i) < 0.02, f(i) = 0.02; end
    % dCP
    fM = f(i);
    dCP(i) = K0 + K1*(1-fM) + K2*sin(pi*fM^kappa);
    if alpha(i) > alpha2_0
                dCP(i) = dCP(i) + (alpha(i)-alpha2_0)*5*K1*exp(3*fM^(-1/kappa));
%         dCP(i) = K0 + K1*exp(-(4*K2/fM)^1);
    end
end

figure; plot(f,0.25-dCP); grid; xlabel('f'); ylabel('$x_{cp}$'); ylim([0.25 0.5])
% figure; plot(alpha*180/pi,dCP); grid; xlabel('alpha, deg'); ylabel('x_{cp}')