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
case_now = 1;
switch case_now
    case 1
        filename{1} = 'gust_test_31.fig';
        filename{2} = 'gust_test_21.fig';
        filename{3} = 'gust_test_22.fig';
        filename{4} = 'gust_test_23.fig';
        filename{5} = 'gust_test_25.fig';
        filename{6} = 'gust_test_27.fig';
        XLim = [0 6];
        lgd_str = {'$\lambda_g$ = -1','$\lambda_g$ = 0','$\lambda_g$ = 1/4','$\lambda_g$ = 1/2','$\lambda_g$ = 1','$\lambda_g$ = 2'};
        authors_ref1 = 'Analytical - Leishman (1997)';
        ylabel = 'Normalized $c_l$';    
end
N = length(filename);
c = lines(N);
markers = {'o','d','s','x','p','*'};
p = gobjects(N,1);

%% Plot compound figure
for i=1:N
    uiopen(filename{i},1); drawnow;
    h = gcf; 
    axes = get(h, 'Children');  
    data = get(axes, 'Children'); data = data(1);
    xdata_mod = data.Children(2).Children(2).XData; ydata_mod = data.Children(2).Children(2).YData;
    xdata_ref1 = data.Children(2).Children(1).XData; ydata_ref1 = data.Children(2).Children(1).YData;
    hold(ax,'on');
    plot(xdata_mod,ydata_mod,'-','Color',c(i,:),'LineWidth',lw,'Parent',ax); 
    plot(xdata_ref1,ydata_ref1,'LineStyle','none','Marker',markers{i},'Color',c(i,:),'MarkerSize',ms,'Parent',ax);
    p(i) = plot(nan,nan,'-','Marker',markers{i},'Color',c(i,:),'LineWidth',lw,'MarkerSize',ms,'Parent',ax);
end
ax.XLim = XLim; ax.XTick = 0:1:6;
ax.FontSize = axes_size; ax.FontName = 'times new roman';
ax.XGrid = 'on'; ax.YGrid = 'on';
ax.XLabel.String = '$\tau$ [semi-chords]'; ax.XLabel.FontWeight = 'normal'; ax.XLabel.FontSize = axes_size;
ax.YLabel.String = ylabel; ax.YLabel.FontWeight = 'normal'; ax.YLabel.FontSize = axes_size;
ax.Title.FontWeight = 'normal'; ax.Title.FontSize = axes_size;
lgd = legend(p,lgd_str); lgd.FontSize = axes_size*0.8;
title(lgd,{'Solid lines: Present model',"Symbols: " + authors_ref1})
