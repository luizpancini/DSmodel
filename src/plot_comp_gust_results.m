clc
clear
close all

% Default plot options
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultLegendInterpreter','latex')
axes_size = 20;
lw = 1;
ms = 5;
% Create figure and tab group
fig = figure('InvertHardcopy','off','Color',[1 1 1]); fig.PaperSize = [8 7.2];
fig.Position = [3.5540e+02 129 920 6.2960e+02];
ax = axes('Parent',fig);
% Select case
case_now = 1;
switch case_now
    case 1
        filename{1} = 'gust_test_1.fig';
        filename{2} = 'gust_test_2.fig';
        filename{3} = 'gust_test_3.fig';
        gust_profile = "1-cos";
        M = 0.2;
        alpha_g = 1.0;
        H = 1*pi;
        XLim = [0 30];
        lgd_str = {'$\alpha_0$ = $0^\circ$','$\alpha_0$ = $10^\circ$','$\alpha_0$ = $15^\circ$'};
        authors_ref = 'CFD - Mallik and Raveh (2019)';
    case 2
        filename{1} = 'gust_test_4.fig';
        filename{2} = 'gust_test_5.fig';
        filename{3} = 'gust_test_6.fig';
        gust_profile = "1-cos";
        M = 0.2;
        alpha_g = 1.0;
        H = 8*pi;
        XLim = [0 80];
        lgd_str = {'$\alpha_0$ = $0^\circ$','$\alpha_0$ = $10^\circ$','$\alpha_0$ = $15^\circ$'};
        authors_ref = 'CFD - Mallik and Raveh (2019)';
    case 3
        filename{1} = 'gust_test_14.fig';
        filename{2} = 'gust_test_18.fig';
        filename{3} = 'gust_test_20.fig';
        gust_profile = "sharp-edge";
        M = 0.5;
        alpha_g = 1.0;
        H = inf;
        XLim = [0 10];
        lgd_str = {'$\lambda_g$ = 0.5','$\lambda_g$ = 1.0','$\lambda_g$ = 2.0'};
        authors_ref = 'CFD - Singh and Baeder (1997)';
end
c = lines(3);
markers = {'o','d','s'};
p = gobjects(3,1);

%% Plot compound figure
for i=1:3
    uiopen(filename{i},1); drawnow;
    h = gcf; 
    axes = get(h, 'Children');  
    data = get(axes, 'Children');
    xdata_mod = get(data{2}(2),'XData'); ydata_mod = get(data{2}(2),'YData');
    xdata_ref = get(data{2}(1),'XData'); ydata_ref = get(data{2}(1),'YData');
    hold(ax,'on');
    plot(xdata_mod,ydata_mod,'-','Color',c(i,:),'LineWidth',lw,'Parent',ax); 
    plot(xdata_ref,ydata_ref,'LineStyle', 'none','Marker',markers{i},'Color',c(i,:),'MarkerSize',ms,'Parent',ax);
    p(i) = plot(nan,nan,'-','Marker',markers{i},'Color',c(i,:),'LineWidth',lw,'MarkerSize',ms,'Parent',ax);
end
ax.XLim = XLim;
ax.FontSize = axes_size; ax.FontName = 'times new roman';
ax.XGrid = 'on'; ax.YGrid = 'on';
ax.XLabel.String = '$\tau$ [semi-chords]'; ax.XLabel.FontWeight = 'normal'; ax.XLabel.FontSize = axes_size;
ax.YLabel.String = '$c_l$ increment'; ax.YLabel.FontWeight = 'normal'; ax.YLabel.FontSize = axes_size;
% ax.Title.String = "Gust profile: " + gust_profile + ", $M$ = " +  num2str(M,'%.2f') + ", $\alpha_g$ = " + num2str(alpha_g,'%.1f') + "$^{\circ}$, $H$ = " + num2str(H,'%.0f'); ax.Title.FontWeight = 'normal'; ax.Title.FontSize = axes_size;
lgd = legend(p,lgd_str); lgd.FontSize = axes_size*0.8;
title(lgd,{'Lines: Model',"Symbols: " + authors_ref})
