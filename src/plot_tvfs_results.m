clc
clear
close all

% Default plot options
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultLegendInterpreter','latex')
axes_size = 20;
lw = 1;
ms = 8;
% Create figure and tab group
fig = figure('InvertHardcopy','off','Color',[1 1 1]); fig.PaperSize = [8 7.2];
fig.Position = [3.5540e+02 129 920 6.2960e+02];
ax = axes('Parent',fig);
% Select case
case_now = 2;
switch case_now
    case 1
        filename{1} = 'tvfs_test_11_cl.fig';
        filename{2} = 'tvfs_test_12_cl.fig';
        filename{3} = 'tvfs_test_13_cl.fig';
        filename{4} = 'tvfs_test_14_cl.fig';
        M = 0.3;
        k = 0.2;
        k_U = 0.2;
        a_0 = 1;
        a_1 = 1;
        psi = 0;
        XLim = [0 360];
        lgd_str = {'$\lambda_U$ = 0.2','$\lambda_U$ = 0.4','$\lambda_U$ = 0.6','$\lambda_U$ = 0.8'};
        authors_ref1 = 'Model - Jose (2006)';
        authors_ref2 = 'CFD - Jose (2006)';
        ylabel = '$c_l / c_{l{qs}}$';
    case 2
        filename{1} = 'tvfs_test_11_cm.fig';
        filename{2} = 'tvfs_test_12_cm.fig';
        filename{3} = 'tvfs_test_13_cm.fig';
        filename{4} = 'tvfs_test_14_cm.fig';
        M = 0.3;
        k = 0.2;
        k_U = 0.2;
        a_0 = 1;
        a_1 = 1;
        psi = 0;
        XLim = [0 360];
        lgd_str = {'$\lambda_U$ = 0.2','$\lambda_U$ = 0.4','$\lambda_U$ = 0.6','$\lambda_U$ = 0.8'};
        authors_ref1 = 'Model - Jose (2006)';
        authors_ref2 = 'CFD - Jose (2006)';
        ylabel = '$c_m$';    
end
c = lines(4);
markers = {'o','d','s','x'};
p = gobjects(4,1);

%% Plot compound figure
for i=1:4
    uiopen(filename{i},1); drawnow;
    h = gcf; 
    axes = get(h, 'Children');  
    data = get(axes, 'Children');
    xdata_mod = get(data{2}(3),'XData'); ydata_mod = get(data{2}(3),'YData');
    xdata_ref1 = get(data{2}(1),'XData'); ydata_ref1 = get(data{2}(1),'YData');
    xdata_ref2 = get(data{2}(2),'XData'); ydata_ref2 = get(data{2}(2),'YData');
    hold(ax,'on');
    plot(xdata_mod,ydata_mod,'-','Color',c(i,:),'LineWidth',lw,'Parent',ax); 
    plot(xdata_ref1,ydata_ref1,'--','Color',c(i,:),'LineWidth',lw,'Parent',ax);
    plot(xdata_ref2,ydata_ref2,'LineStyle', 'none','Marker',markers{i},'Color',c(i,:),'MarkerSize',ms,'Parent',ax);
    p(i) = plot(nan,nan,'-','Marker',markers{i},'Color',c(i,:),'LineWidth',lw,'MarkerSize',ms,'Parent',ax);
end
ax.XLim = XLim; ax.XTick = 0:45:360;
ax.FontSize = axes_size; ax.FontName = 'times new roman';
ax.XGrid = 'on'; ax.YGrid = 'on';
ax.XLabel.String = '$\omega t$ [deg]'; ax.XLabel.FontWeight = 'normal'; ax.XLabel.FontSize = axes_size;
ax.YLabel.String = ylabel; ax.YLabel.FontWeight = 'normal'; ax.YLabel.FontSize = axes_size;
ax.Title.String = "$M_0$ = " +  num2str(M,'%.2f') + ", $k_U$ = " + num2str(k_U,'%.2f') + ", $\alpha_0$ = " +  num2str(a_0,'%.1f') + "$^{\circ}$, $\alpha_1$ = " + num2str(a_1,'%.1f') + "$^{\circ}$" + ", $k$ = " + num2str(k,'%.2f') + ", $\psi$ = " + num2str(psi,'%.0f') + "$^{\circ}$";
ax.Title.FontWeight = 'normal'; ax.Title.FontSize = axes_size;
lgd = legend(p,lgd_str); lgd.FontSize = axes_size*0.8;
title(lgd,{'Solid lines: Present model',"Broken lines: " + authors_ref1,"Symbols: " + authors_ref2})
