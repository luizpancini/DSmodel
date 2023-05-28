clc
clear
close all

% Constants
f0 = 0.035;
fb1 = 0.375;
fb2 = 0.007;
alpha1 = 11.25*pi/180;
alpha2 = 21.0*pi/180;
S1 = 4.0*pi/180;
S2 = 4.0*pi/180;
S3 = 1.3;

% Load experimental data
filepath = sprintf('../OSU Data/S809/S809_f_vs_alpha_Re_1e6.mat');  
load(filepath)
alpha_exp = S809_f_vs_alpha_Re_1e6(:,1);
f_exp = S809_f_vs_alpha_Re_1e6(:,2);
clear S809_f_vs_alpha_Re_1e6

% Set AoA vector
alpha = pi/180*(0:0.01:40);

% Get separation point
f = nan(length(alpha),1);
for i=1:length(alpha)
    if alpha(i) <= alpha1
        f_i = 1-(1-fb1)*exp((alpha(i)-alpha1)/S1);
    elseif alpha(i) > alpha1 && alpha(i) < alpha2
        f_i = fb2+(fb1-fb2)*((alpha2-alpha(i))/(alpha2-alpha1))^S3;
    else
        f_i = f0+(fb2-f0)*exp((alpha2-alpha(i))/S2);
    end
    f(i) = f_i;
end

% Plot
figure;
plot(alpha*180/pi,f,'k-',alpha_exp,f_exp,'ko'); ylim([0,1]);
grid