%% Get data
% Dynamic
alpha_cl_exp =  OUTPUTS.data.alpha_exp_cl;
alpha_cm_exp =  OUTPUTS.data.alpha_exp_cm;
alpha_cd_exp =  OUTPUTS.data.alpha_exp_cd;
time_exp_cl =   OUTPUTS.data.time_exp_cl; 
time_exp_cm =   OUTPUTS.data.time_exp_cm; 
time_exp_cd =   OUTPUTS.data.time_exp_cd; 
clt_exp =       OUTPUTS.data.clt_exp;
cmt_exp =       OUTPUTS.data.cmt_exp;
cdt_exp =       OUTPUTS.data.cdt_exp;
cl_exp =        OUTPUTS.data.cl_exp;
cm_exp =        OUTPUTS.data.cm_exp;
cd_exp =        OUTPUTS.data.cd_exp;
alpha_mod =     OUTPUTS.outputs.alpha;
cl_mod =        OUTPUTS.outputs.c_l;
cm_mod =        OUTPUTS.outputs.c_m;
cd_mod =        OUTPUTS.outputs.c_d;

time_exp_cl = (time_exp_cl+90)/360; i=find(time_exp_cl>=0.25,1,'first'); time_exp_cl(1:i-1) = time_exp_cl(1:i-1)+1; time_exp_cl = circshift(time_exp_cl,-(i-1)); time_exp_cl = time_exp_cl - time_exp_cl(1) + time_exp_cl(end-i+2)-1;
clt_exp = circshift(clt_exp,-(i-1));
time_exp_cm = (time_exp_cm+90)/360; i=find(time_exp_cm>=0.25,1,'first'); time_exp_cm(1:i-1) = time_exp_cm(1:i-1)+1; time_exp_cm = circshift(time_exp_cm,-(i-1)); time_exp_cm = time_exp_cm - time_exp_cm(1) + time_exp_cm(end-i+2)-1;
cmt_exp = circshift(cmt_exp,-(i-1));
time_exp_cd = (time_exp_cd+90)/360; i=find(time_exp_cd>=0.25,1,'first'); time_exp_cd(1:i-1) = time_exp_cd(1:i-1)+1; time_exp_cd = circshift(time_exp_cd,-(i-1)); time_exp_cd = time_exp_cd - time_exp_cd(1) + time_exp_cd(end-i+2)-1;
cdt_exp = circshift(cdt_exp,-(i-1));

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
c_mI = OUTPUTS.outputs.c_mI;
t_cycle = OUTPUTS.params.t_cycle; itt = OUTPUTS.interp.itt; ittf = OUTPUTS.interp.ittf;
t_plot = (OUTPUTS.outputs.t(itt:ittf)-OUTPUTS.outputs.t(itt))/t_cycle;

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
alpha_cl_static_exp = h_cl.CurrentAxes.Children(1).XData;
cl_static_exp = h_cl.CurrentAxes.Children(1).YData;
cl_static = h_cl.CurrentAxes.Children(2).YData;
close all
uiopen('cm_static_BL.fig',1)
h_cm = gcf; 
alpha_cm_static_exp = h_cm.CurrentAxes.Children(1).XData;
cm_static_exp = h_cm.CurrentAxes.Children(1).YData;
cm_static = h_cm.CurrentAxes.Children(2).YData;
close all
uiopen('cd_static_BL.fig',1)
h_cd = gcf; 
alpha_cd_static_exp = h_cd.CurrentAxes.Children(1).XData;
cd_static_exp = h_cd.CurrentAxes.Children(1).YData;
alpha_static = h_cd.CurrentAxes.Children(2).XData;
cd_static = h_cd.CurrentAxes.Children(2).YData;
close all

% Unsteady - UM/NAST static stall
alpha_stall = 13.7*pi/180;
% cl_stall = 1.365;
for i=1:length(alpha_C)
    cl_unsteady(i) = c_n_alpha/sqrt(1-0.3^2).*min([alpha_C(i),alpha_stall])+c_nI(i);
end
c_m_alpha = 0.035/(13.6*pi/180);
cm_stall = 0.029;
cm_unsteady = c_m0+alpha*c_m_alpha;
cm_unsteady(cm_unsteady>cm_stall) = cm_stall;
cm_unsteady = cm_unsteady+c_mI*1.7;
figure;plot(alpha*180/pi,cl_unsteady);grid

%% Generate cl animated plot
clc; close all
% Default plot options
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultLegendInterpreter','latex')
axes_size = 20;
lw = 1.5;
ms = 6;
% Create figure and ajust settings
fig = figure('InvertHardcopy','off','Color',[1 1 1]); fig.PaperSize = [9 7.2];
fig.Position = [1.5860e+02   9.7000e+01   1.2072e+03   7.4400e+02];
ax = axes('Parent',fig); hold(ax,'on');
ax.XLim = [0, 25]; ax.YLim = [0, 2.0];
ax.FontSize = axes_size; ax.FontName = 'times new roman';
ax.XGrid = 'on'; ax.YGrid = 'on';
ax.XLabel.String = 'Angle of attack [deg]'; ax.XLabel.FontWeight = 'normal'; ax.XLabel.FontSize = axes_size;
ax.YLabel.String = 'Lift coefficient'; ax.YLabel.FontWeight = 'normal'; ax.YLabel.FontSize = axes_size;
pause_time = 2e-1;
save_gif = 1;
gif_filename = 'cl_vs_alpha.gif';
[Y,Z] = import_airfoil_points('NACA0012.dat');
Y = Y-(0.25);
p1 = gobjects(4,1);
N = round(length(alpha)/4)+1;
% Plot curves
plot(alpha_cl_static_exp,cl_static_exp,'ko','LineWidth',lw,'MarkerSize',ms,'Parent',ax);
plot(alpha_static,cl_static,'k-','LineWidth',lw,'MarkerSize',ms,'Parent',ax);
p = plot(nan,nan,'ko',nan,nan,'k-','LineWidth',lw,'MarkerSize',ms,'Parent',ax);
lgd = legend(p,{'Static - Exp','Static - DSM'}); lgd.AutoUpdate = 'off'; lgd.FontSize = axes_size*0.8; lgd.Position = [6.9574e-01   1.2622e-01   2.3353e-01   2.1283e-01];
plot(alpha_cl_exp,cl_exp,'bo','LineWidth',lw,'MarkerSize',ms,'Parent',ax);
p = plot(nan,nan,'ko',nan,nan,'k-',nan,nan,'bo','LineWidth',lw,'MarkerSize',ms,'Parent',ax);
lgd = legend(p,{'Static - Exp','Static - DSM','Dynamic - Exp'});
ax.Title.String = "$\alpha = 12.0^\circ + 10.0^\circ$ sin$(kUt/b), k \approx 0.1, U \approx 100 m/s$";
% pause(5);
h2 = animatedline(ax,alpha(1)*180/pi,cl_mod(1),'Color','b','LineWidth',lw); 
p = plot(nan,nan,'ko',nan,nan,'k-',nan,nan,'bo',nan,nan,'b-','LineWidth',lw,'MarkerSize',ms,'Parent',ax);
lgd = legend(p,{'Static - Exp','Static - DSM','Dynamic - Exp','Dynamic - DSM'}); 
ax2 = axes('Position',[1.6200e-01   7.5000e-01   3.3539e-01   1.7503e-01],'Box','on'); hold(ax2,'on');
ax2.XLim = [0, 1]; ax2.YLim = ax.YLim; ax2.XGrid = 'on'; ax2.YGrid = 'on'; ax2.XTick = [0, 0.25 0.5 0.75 1]; ax2.XTickLabel = {'0','1/4','1/2','3/4','1'};
h2t = animatedline(ax2,t_plot(1),cl_mod(1),'Color','b','LineWidth',lw);
ax3 = axes('Position',[3.9243e-01   1.2865e-01   2.5017e-01   1.1150e-01],'Box','on'); hold(ax3,'off'); 
ax2.Visible = 'off';
if save_gif, drawnow; gif(gif_filename,'overwrite',true); end
for k = 1:N
    delete(p1);
    % Plot in AoA domain
    addpoints(h2,alpha(k)*180/pi,cl_mod(k));
    p1(1) = plot(alpha(k)*180/pi,cl_mod(k),'b*','LineStyle','none','MarkerSize',2*ms,'Parent',ax);
    % Plot in time domain
%     addpoints(h2t,t_plot(k),cl_mod(k));
%     p1(2) = plot(t_plot(k),cl_mod(k),'b*','LineStyle','none','MarkerSize',ms,'Parent',ax2);
    % Plot airfoil
    R = R_E321([0,0,alpha(k)]);
    rot_coords = R*[zeros(1,length(Y)); Y'; Z'];
    Yr = rot_coords(2,:)';
    Zr = rot_coords(3,:)';
    plot(Yr+.25,-Zr,'k-',0.25,0,'ko','LineWidth',1,'MarkerSize',6,'Parent',ax3);
    ax3.XLim = [-0.2 1.2]; ax3.YLim = [-0.30 0.20];
    ax3.XTickLabel = {}; ax3.YTickLabel = {}; ax3.DataAspectRatio = [1 1 1]; 
    ax3.Title.String = "$\alpha = $" + num2str(180/pi*alpha(k),'%.1f') + "$^ \circ$"; ax3.Title.FontSize = 16; ax3.Title.Position = [8.8712e-01   1.1249e-03  0];
    % Set textarrows
    if k == round(N/32)
        annotation(fig,'textarrow',[0.623478260869565 0.568695652173913],[0.71156289707751 0.705209656925032],'String',{'Stall delayed'},'Color','b','FontSize',14);
    elseif k == round(N/16)
        annotation(fig,'textarrow',[0.587826086956522 0.62],[0.851604828462516 0.809402795425667],'String',{'DS onset'},'Color','b','FontSize',14);
    elseif k == round(N/8)
        annotation(fig,'textarrow',[0.784347826086958 0.743913043478265],[0.909783989834816 0.886277001270648],'String',{'DSV peak load'},'Color','b','FontSize',14);
    elseif k == round(N/4)
        annotation(fig,'textarrow',[0.820869565217391 0.809130434782611],[0.625158831003812 0.574968233799238],'String',{'2nd DSV'},'Color','b','FontSize',14);
    elseif k == round(N/4*1.1)
        annotation(fig,'textarrow',[0.836521739130435 0.824347826086956],[0.467598475222363 0.50571791613723],'String',{'Stalled flow','convection'},'Color','b','FontSize',14);
    elseif k == round(N/2*1.1)
        annotation(fig,'textarrow',[0.351304347826087 0.378260869565217],[0.231257941550191 0.29987293519695],'String',{'Flow reattachment','begins'},'Color','b','FontSize',14);
    elseif k == round(N*0.75)
        annotation(fig,'textarrow',[0.185217391304348 0.188695652173913],[0.348157560355781 0.238881829733164],'String',{'Attached flow'},'Color','b','FontSize',14);
    end
    % Draw and save gif
    drawnow 
    if save_gif
        gif('overwrite',true,'DelayTime',pause_time);
    end
end
gif_filename = 'cl_vs_alpha2.gif';
if save_gif, drawnow; gif(gif_filename,'overwrite',true); end
h2 = animatedline(ax,alpha(1)*180/pi,cl_mod(1),'Color','b','LineWidth',lw);
h3 = animatedline(ax,alpha(1)*180/pi,cl_unsteady(1),'Color','r','LineWidth',lw);
h2t = animatedline(ax2,t_plot(1),cl_mod(1),'Color','b','LineWidth',lw);
h3t = animatedline(ax2,t_plot(1),cl_unsteady(1),'Color','r','LineWidth',lw);
p = plot(nan,nan,'ko',nan,nan,'k-',nan,nan,'bo',nan,nan,'b-',nan,nan,'r-','LineWidth',lw,'MarkerSize',ms,'Parent',ax);
lgd = legend(p,{'Static - Exp','Static - DSM','Dynamic - Exp','Dynamic - DSM','Dynamic - NAST'});
plot(time_exp_cl,clt_exp,'bo','LineWidth',lw,'MarkerSize',ms/2,'Parent',ax2); ax2.Visible = 'on';
for k = 1:N
    delete(p1);
    % Plot in AoA domain
    addpoints(h2,alpha(k)*180/pi,cl_mod(k));
    addpoints(h3,alpha(k)*180/pi,cl_unsteady(k));
    p1(1) = plot(alpha(k)*180/pi,cl_mod(k),'b*','LineStyle','none','MarkerSize',2*ms,'Parent',ax);
    p1(2) = plot(alpha(k)*180/pi,cl_unsteady(k),'r*','LineStyle','none','MarkerSize',2*ms,'Parent',ax);
    % Plot in time domain
    addpoints(h2t,t_plot(k),cl_mod(k));
    addpoints(h3t,t_plot(k),cl_unsteady(k));
    p1(3) = plot(t_plot(k),cl_mod(k),'b*','LineStyle','none','MarkerSize',ms,'Parent',ax2);
    p1(4) = plot(t_plot(k),cl_unsteady(k),'r*','LineStyle','none','MarkerSize',ms,'Parent',ax2);
    % Plot airfoil
    R = R_E321([0,0,alpha(k)]);
    rot_coords = R*[zeros(1,length(Y)); Y'; Z'];
    Yr = rot_coords(2,:)';
    Zr = rot_coords(3,:)';
    plot(Yr+.25,-Zr,'k-',0.25,0,'ko','LineWidth',1,'MarkerSize',6,'Parent',ax3);
    ax3.XLim = [-0.2 1.2]; ax3.YLim = [-0.30 0.20];
    ax3.XTickLabel = {}; ax3.YTickLabel = {}; ax3.DataAspectRatio = [1 1 1]; 
    ax3.Title.String = "$\alpha = $" + num2str(180/pi*alpha(k),'%.1f') + "$^ \circ$"; ax3.Title.FontSize = 16; ax3.Title.Position = [8.8712e-01   1.1249e-03  0];
    % Draw and save gif
    drawnow 
    if save_gif
        gif('overwrite',true,'DelayTime',pause_time/2);
    end
end

%% Generate cm animated plot
% Create figure and ajust settings
fig = figure('InvertHardcopy','off','Color',[1 1 1]); fig.PaperSize = [9 7.2];
fig.Position = [1.5860e+02   9.7000e+01   1.2072e+03   7.4400e+02];
ax = axes('Parent',fig); hold(ax,'on');
ax.XLim = [0, 25]; ax.YLim = [-0.4, 0.1];
ax.FontSize = axes_size; ax.FontName = 'times new roman';
ax.XGrid = 'on'; ax.YGrid = 'on';
ax.XLabel.String = 'Angle of attack [deg]'; ax.XLabel.FontWeight = 'normal'; ax.XLabel.FontSize = axes_size;
ax.YLabel.String = 'Moment coefficient'; ax.YLabel.FontWeight = 'normal'; ax.YLabel.FontSize = axes_size;
pause_time = 1e-1;
save_gif = 1;
gif_filename = 'cm_vs_alpha.gif';
[Y,Z] = import_airfoil_points('NACA0012.dat');
Y = Y-(0.25);
p1 = gobjects(4,1);
N = round(length(alpha)/4)+1;
% Plot curves
plot(alpha_cm_static_exp,cm_static_exp,'ko','LineWidth',lw,'MarkerSize',ms,'Parent',ax);
plot(alpha_static,cm_static,'k-','LineWidth',lw,'MarkerSize',ms,'Parent',ax);
p = plot(nan,nan,'ko',nan,nan,'k-','LineWidth',lw,'MarkerSize',ms,'Parent',ax);
lgd = legend(p,{'Static - Exp','Static - DSM'}); lgd.AutoUpdate = 'off'; lgd.FontSize = axes_size*0.8; lgd.Position = [1.5896e-01   4.0364e-01   2.3353e-01   2.1283e-01];
plot(alpha_cm_exp,cm_exp,'bo','LineWidth',lw,'MarkerSize',ms,'Parent',ax);
p = plot(nan,nan,'ko',nan,nan,'k-',nan,nan,'bo','LineWidth',lw,'MarkerSize',ms,'Parent',ax);
lgd = legend(p,{'Static - Exp','Static - DSM','Dynamic - Exp'});
ax.Title.String = "$\alpha = 12.0^\circ + 10.0^\circ$ sin$(kUt/b), k \approx 0.1, U \approx 100 m/s$";
% pause(5);
h2 = animatedline(ax,alpha(1)*180/pi,cm_mod(1),'Color','b','LineWidth',lw); 
p = plot(nan,nan,'ko',nan,nan,'k-',nan,nan,'bo',nan,nan,'b-','LineWidth',lw,'MarkerSize',ms,'Parent',ax);
lgd = legend(p,{'Static - Exp','Static - DSM','Dynamic - Exp','Dynamic - DSM'}); 
ax2 = axes('Position',[1.6465e-01   1.4785e-01   3.3539e-01   1.7503e-01],'Box','on'); hold(ax2,'on');
ax2.XLim = [0, 1]; ax2.YLim = ax.YLim; ax2.XGrid = 'on'; ax2.YGrid = 'on'; ax2.XTick = [0, 0.25 0.5 0.75 1]; ax2.XTickLabel = {'0','1/4','1/2','3/4','1'};
h2t = animatedline(ax2,t_plot(1),cm_mod(1),'Color','b','LineWidth',lw);
ax2.Visible = 'off';
ax3 = axes('Position',[4.2557e-01   3.7919e-01   2.5017e-01   1.1150e-01],'Box','on'); hold(ax3,'off'); 
if save_gif, drawnow; gif(gif_filename,'overwrite',true); end
for k = 1:N
    delete(p1);
    % Plot in AoA domain
    addpoints(h2,alpha(k)*180/pi,cm_mod(k));
    p1(1) = plot(alpha(k)*180/pi,cm_mod(k),'b*','LineStyle','none','MarkerSize',2*ms,'Parent',ax);
    % Plot in time domain
%     addpoints(h2t,t_plot(k),cm_mod(k));
%     p1(2) = plot(t_plot(k),cm_mod(k),'b*','LineStyle','none','MarkerSize',ms,'Parent',ax2);
    % Plot airfoil
    R = R_E321([0,0,alpha(k)]);
    rot_coords = R*[zeros(1,length(Y)); Y'; Z'];
    Yr = rot_coords(2,:)';
    Zr = rot_coords(3,:)';
    plot(Yr+.25,-Zr,'k-',0.25,0,'ko','LineWidth',1,'MarkerSize',6,'Parent',ax3);
    ax3.XLim = [-0.2 1.2]; ax3.YLim = [-0.30 0.20];
    ax3.XTickLabel = {}; ax3.YTickLabel = {}; ax3.DataAspectRatio = [1 1 1]; 
    ax3.Title.String = "$\alpha = $" + num2str(180/pi*alpha(k),'%.1f') + "$^ \circ$"; ax3.Title.FontSize = 16; ax3.Title.Position = [8.8712e-01   1.1249e-03  0];
    % Set textarrows
    if k == round(N/32)
        annotation(fig,'textarrow',[0.615639496355202 0.60569913850232],[0.598924731182796 0.70215053763441],'String',{'Stall delayed'},'Color','b','FontSize',14);
    elseif k == round(N/16)
        annotation(fig,'textarrow',[0.684559310801856 0.664015904572565],[0.783870967741936 0.73225806451613],'String',{'DS onset'},'Color','b','FontSize',14);
    elseif k == round(N/8)
        annotation(fig,'textarrow',[0.754141815772034 0.746189529489728],[0.195698924731183 0.258064516129032],'String',{'DSV peak load'},'Color','b','FontSize',14);
    elseif k == round(N/4)
        annotation(fig,'textarrow',[0.825049701789264 0.809807819748178],[0.437634408602151 0.493548387096775],'String',{'2nd DSV'},'Color','b','FontSize',14);
    elseif k == round(N/4*1.1)
        annotation(fig,'textarrow',[0.833664678595096 0.811795891318754],[0.617204301075269 0.54516129032258],'String',{'Stalled flow','convection'},'Color','b','FontSize',14);
    elseif k == round(N/2*1.1)
        annotation(fig,'textarrow',[0.346587143803843 0.408880053015241],[0.891397849462368 0.869892473118283],'String',{'Flow reattachment','begins'},'Color','b','FontSize',14);
    elseif k == round(N*0.75)
        annotation(fig,'textarrow',[0.204108681245858 0.201457919151756],[0.648387096774194 0.730107526881721],'String',{'Attached flow'},'Color','b','FontSize',14);
    end
    % Draw and save gif
    drawnow 
    pause(pause_time);
    if save_gif
        gif('overwrite',true);
    end
end
h2 = animatedline(ax,alpha(1)*180/pi,cm_mod(1),'Color','b','LineWidth',lw);
h3 = animatedline(ax,alpha(1)*180/pi,cm_unsteady(1),'Color','r','LineWidth',lw);
h2t = animatedline(ax2,t_plot(1),cm_mod(1),'Color','b','LineWidth',lw);
h3t = animatedline(ax2,t_plot(1),cm_unsteady(1),'Color','r','LineWidth',lw);
p = plot(nan,nan,'ko',nan,nan,'k-',nan,nan,'bo',nan,nan,'b-',nan,nan,'r-','LineWidth',lw,'MarkerSize',ms,'Parent',ax);
lgd = legend(p,{'Static - Exp','Static - DSM','Dynamic - Exp','Dynamic - DSM','Dynamic - NAST'});
plot(time_exp_cm,cmt_exp,'bo','LineWidth',lw,'MarkerSize',ms/2,'Parent',ax2);
ax2.Visible = 'on';
for k = 1:N
    delete(p1);
    % Plot in AoA domain
    addpoints(h2,alpha(k)*180/pi,cm_mod(k));
    addpoints(h3,alpha(k)*180/pi,cm_unsteady(k));
    p1(1) = plot(alpha(k)*180/pi,cm_mod(k),'b*','LineStyle','none','MarkerSize',2*ms,'Parent',ax);
    p1(2) = plot(alpha(k)*180/pi,cm_unsteady(k),'r*','LineStyle','none','MarkerSize',2*ms,'Parent',ax);
    % Plot in time domain
    addpoints(h2t,t_plot(k),cm_mod(k));
    addpoints(h3t,t_plot(k),cm_unsteady(k));
    p1(3) = plot(t_plot(k),cm_mod(k),'b*','LineStyle','none','MarkerSize',ms,'Parent',ax2);
    p1(4) = plot(t_plot(k),cm_unsteady(k),'r*','LineStyle','none','MarkerSize',ms,'Parent',ax2);
    % Plot airfoil
    R = R_E321([0,0,alpha(k)]);
    rot_coords = R*[zeros(1,length(Y)); Y'; Z'];
    Yr = rot_coords(2,:)';
    Zr = rot_coords(3,:)';
    plot(Yr+.25,-Zr,'k-',0.25,0,'ko','LineWidth',1,'MarkerSize',6,'Parent',ax3);
    ax3.XLim = [-0.2 1.2]; ax3.YLim = [-0.30 0.20];
    ax3.XTickLabel = {}; ax3.YTickLabel = {}; ax3.DataAspectRatio = [1 1 1]; 
    ax3.Title.String = "$\alpha = $" + num2str(180/pi*alpha(k),'%.1f') + "$^ \circ$"; ax3.Title.FontSize = 16; ax3.Title.Position = [8.8712e-01   1.1249e-03  0];
    % Draw and save gif
    drawnow 
    pause(pause_time);
    if save_gif
        gif('overwrite',true);
    end
end