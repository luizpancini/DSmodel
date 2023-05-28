%% Get data
% Dynamic
alpha_cl_exp =  OUTPUTS.data.alpha_exp_cl;
alpha_cm_exp =  OUTPUTS.data.alpha_exp_cm;
alpha_cd_exp =  OUTPUTS.data.alpha_exp_cd;
cl_exp =        OUTPUTS.data.cl_exp;
cm_exp =        OUTPUTS.data.cm_exp;
cd_exp =        OUTPUTS.data.cd_exp;
alpha_mod =     OUTPUTS.outputs.alpha;
cl_mod =        OUTPUTS.outputs.c_l;
cm_mod =        OUTPUTS.outputs.c_m;
cd_mod =        OUTPUTS.outputs.c_d;

% Static - model
alpha_0L = OUTPUTS.params.alpha_0L;
alpha1_0c = OUTPUTS.params.alpha1_0c;
kappa_0 = OUTPUTS.params.kappa_0;
eta = OUTPUTS.params.eta;
c_d0 = OUTPUTS.params.c_d0;
c_m0 = OUTPUTS.params.c_m0;
c_n_alpha = OUTPUTS.params.c_n_alpha;
K0 = OUTPUTS.params.K0;
K1 = OUTPUTS.params.K1;
K2 = OUTPUTS.params.K2;
U = OUTPUTS.outputs.U;
w_QS = OUTPUTS.outputs.w_QS;
alpha = OUTPUTS.outputs.alpha;
alpha_C = OUTPUTS.outputs.alpha_C;
f_n = OUTPUTS.outputs.f_n;
f_m = OUTPUTS.outputs.f_m;
f_c = OUTPUTS.outputs.f_c;
c_nI = OUTPUTS.outputs.c_nI;
cn_static_mod = c_n_alpha./U.*(w_QS-U.*sin(alpha_0L)).*((1+sqrt(f_n))/2).^2;
c_mf = cn_static_mod.*(K0+K1*(1-f_m)+K2*sin(pi*f_m.^kappa_0));
cm_static_mod = c_m0+c_mf;
cc_static_mod = -c_d0*cos(alpha)+c_n_alpha*sin(alpha_C-alpha_0L).^2.*f_c.^(1/2+abs(alpha)./alpha1_0c);
R = zeros(2,2,length(alpha));
for i=1:length(alpha)
    R(:,:,i) = [cos(alpha(i)),sin(alpha(i));sin(alpha(i)),-cos(alpha(i))];
    coefs = R(:,:,i)^-1*[cn_static_mod(i); cc_static_mod(i)];
    cl_static_mod(i) = coefs(1);
    cd_static_mod(i) = coefs(2);
end

% Static - experimental
uiopen('cl_static_BL.fig',1)
h_cl = gcf; 
alpha_cl_static_exp = h_cl.CurrentAxes.Children(2).XData;
cl_static_exp = h_cl.CurrentAxes.Children(2).YData;
uiopen('cm_static_BL.fig',1)
h_cm = gcf; 
alpha_cm_static_exp = h_cm.CurrentAxes.Children(2).XData;
cm_static_exp = h_cm.CurrentAxes.Children(2).YData;
uiopen('cd_static_BL.fig',1)
h_cd = gcf; 
alpha_cd_static_exp = h_cd.CurrentAxes.Children(2).XData;
cd_static_exp = h_cd.CurrentAxes.Children(2).YData;
close all

% Unsteady - UM/NAST static stall
cl_unsteady = (c_n_alpha/sqrt(1-0.3^2).*alpha_C+c_nI*0.4); % 0.4 is a factor in this case between incompressible and compressible flow inertial cn
cl_unsteady(cl_unsteady>1.4) = 1.4;

%% Generate animated plot
clc; close all
% Default plot options
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultLegendInterpreter','latex')
axes_size = 20;
lw = 1.5;
ms = 6;
% Create figure and ajust settings
fig = figure('InvertHardcopy','off','Color',[1 1 1]); fig.PaperSize = [9 7.2];
fig.Position = [3.5540e+02 129 920 6.2960e+02];
ax = axes('Parent',fig); hold(ax,'on');
ax.XLim = [0, 25]; ax.YLim = [0, 2.0];
ax.FontSize = axes_size; ax.FontName = 'times new roman';
ax.XGrid = 'on'; ax.YGrid = 'on';
ax.XLabel.String = 'Angle of attack [deg]'; ax.XLabel.FontWeight = 'normal'; ax.XLabel.FontSize = axes_size;
ax.YLabel.String = 'Lift coefficient'; ax.YLabel.FontWeight = 'normal'; ax.YLabel.FontSize = axes_size;
ax.Title.String = "$\alpha = 12.0^\circ + 10.0^\circ$ sin$(kUt/b), k \approx 0.1, U \approx 100 m/s$";
p = plot(nan,nan,'ko',nan,nan,'k-',nan,nan,'bo',nan,nan,'b-',nan,nan,'r-','LineWidth',lw,'MarkerSize',ms,'Parent',ax);
lgd = legend(p,{'Static - Exp','Static - DSM','Dynamic - Exp','Dynamic - DSM','Dynamic - NAST'}); lgd.FontSize = axes_size*0.8; lgd.AutoUpdate = 'off'; lgd.Location = 'southeast';
ax2 = axes('Position',[1.6200e-01   7.5000e-01   3.3539e-01   1.7503e-01],'Box','on'); hold(ax2,'off');
pause_time = 10e-2;
save_gif = 1;
gif_filename = 'cl_vs_alpha.gif';
[Y,Z] = import_airfoil_points('NACA0012.dat');
Y = Y-(0.25);
p1 = gobjects(3,1);
N = round(length(alpha)/4);
% Plot curves
plot(alpha_cl_static_exp(1:200:end),cl_static_exp(1:200:end),'ko','LineWidth',lw,'MarkerSize',ms,'Parent',ax);
plot(alpha_cl_exp,cl_exp,'bo','LineWidth',lw,'MarkerSize',ms,'Parent',ax);
h1 = animatedline(ax,alpha(1)*180/pi,cl_static_mod(1),'Color','k','LineWidth',lw); 
h2 = animatedline(ax,alpha(1)*180/pi,cl_mod(1),'Color','b','LineWidth',lw); 
h3 = animatedline(ax,alpha(1)*180/pi,cl_unsteady(1),'Color','r','LineWidth',lw); 
for k = 1:N
    delete(p1);
    addpoints(h1,alpha(k)*180/pi,cl_static_mod(k));
    addpoints(h2,alpha(k)*180/pi,cl_mod(k));
    addpoints(h3,alpha(k)*180/pi,cl_unsteady(k));
    p1(1) = plot(alpha(k)*180/pi,cl_static_mod(k),'k*','LineStyle','none','MarkerSize',2*ms,'Parent',ax);
    p1(2) = plot(alpha(k)*180/pi,cl_mod(k),'b*','LineStyle','none','MarkerSize',2*ms,'Parent',ax);
    p1(3) = plot(alpha(k)*180/pi,cl_unsteady(k),'r*','LineStyle','none','MarkerSize',2*ms,'Parent',ax);
    R = R_E321([0,0,alpha(k)]);
    rot_coords = R*[zeros(1,length(Y)); Y'; Z'];
    Yr = rot_coords(2,:)';
    Zr = rot_coords(3,:)';
    plot(Yr+.25,-Zr,'k-',0.25,0,'ko','LineWidth',1,'MarkerSize',6,'Parent',ax2);
    ax2.XTickLabel = {}; ax2.YTickLabel = {}; ax2.XGrid = 'off'; ax2.YGrid = 'off'; 
    ax2.DataAspectRatio = [1 1 1]; ax2.Title.String = "$\alpha = $" + num2str(180/pi*alpha(k),'%.1f') + "$^ \circ$"; ax2.Title.FontSize = 16; ax2.Title.Position = [9.0985e-01   2.9534e-02  0];
    ax2.XLim = [-0.2 1.2]; ax2.YLim = [-0.30 0.20];
    if k == round(N/32)
        annotation(fig,'textarrow',[0.623478260869565 0.568695652173913],[0.71156289707751 0.705209656925032],'String',{'Stall delayed'},'Color','b');
    elseif k == round(N/16)
        annotation(fig,'textarrow',[0.587826086956522 0.62],[0.851604828462516 0.809402795425667],'String',{'DS onset'},'Color','b');
    elseif k == round(N/8)
        annotation(fig,'textarrow',[0.784347826086958 0.743913043478265],[0.909783989834816 0.886277001270648],'String',{'DSV peak load'},'Color','b');
    elseif k == round(N/4)
        annotation(fig,'textarrow',[0.820869565217391 0.809130434782611],[0.625158831003812 0.574968233799238],'String',{'2nd DSV'},'Color','b');
    elseif k == round(N/4*1.1)
        annotation(fig,'textarrow',[0.836521739130435 0.824347826086956],[0.467598475222363 0.50571791613723],'String',{'Stalled flow','convection'},'Color','b');
    elseif k == round(N/2*1.1)
        annotation(fig,'textarrow',[0.354782608695652 0.378260869565217],[0.219822109275731 0.29987293519695],'String',{'Flow reattachment begins'},'Color','b');
    elseif k == round(N*0.75)
        annotation(fig,'textarrow',[0.185217391304348 0.188695652173913],[0.348157560355781 0.238881829733164],'String',{'Attached flow'},'Color','b');
    end
    % Draw and save gif
    drawnow 
    pause(pause_time);
    if save_gif
        if k == 1
            gif(gif_filename,'overwrite',true);
        else
            gif('overwrite',true);
        end
    end
end
