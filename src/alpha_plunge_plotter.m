function [h,fig] = alpha_plunge_plotter(tabbed,fig,tab_name,coef_name,alpha_plunge,coef,alpha_plunge_exp,coef_exp,alpha_plunge_mod,coef_mod,authors,iti,tabgp,model_name,plot_opt)

% Check input
if diff(plot_opt.xlim_alpha_plunge_vec) == 0
    h = gobjects(0); return;
end

%% Initialize tab and set axes design
if tabbed
    tab = uitab('Parent',tabgp,'Title',tab_name,'BackgroundColor',[1 1 1]);
    ax = axes('Parent',tab,'FontSize',plot_opt.axes_size,'FontName','times new roman');
    ax.Title.String = plot_opt.title_str; ax.Title.FontWeight = 'normal'; ax.Title.FontSize = plot_opt.axes_size*0.6;
else
    ax = axes('Parent',fig,'FontSize',plot_opt.axes_size,'FontName','times new roman');
end
hold(ax,'on'); ax.XGrid = 'on'; ax.YGrid = 'on';
ax.XLim = plot_opt.xlim_alpha_plunge_vec;
ax.XLabel.String = '$\alpha_{\textnormal{eq}}$ [deg]'; ax.XLabel.FontWeight = 'normal'; ax.XLabel.FontSize = plot_opt.axes_size;
ax.YLabel.String = coef_name; ax.YLabel.FontWeight = 'normal'; ax.YLabel.FontSize = plot_opt.axes_size;

%% Plot
alpha_plunge_plot = alpha_plunge(iti:end)*180/pi; coef_plot = coef(iti:end);
if isnan(alpha_plunge_exp)
    h = plot(alpha_plunge_plot,coef_plot,'k-','LineWidth',plot_opt.lw,'MarkerSize',plot_opt.ms,'Parent',ax);
elseif length(authors) == 1
    h = plot(alpha_plunge_plot,coef_plot,'k-',alpha_plunge_exp,coef_exp,'ko','LineWidth',plot_opt.lw,'MarkerSize',plot_opt.ms,'Parent',ax);
elseif length(authors) == 2
    h = plot(alpha_plunge_plot,coef_plot,'k-',alpha_plunge_exp,coef_exp,'ko',alpha_plunge_mod,coef_mod,'k--','LineWidth',plot_opt.lw,'MarkerSize',plot_opt.ms,'Parent',ax);
end
% Set legend
if isnan(alpha_plunge_exp)
    lgd = legend(h,{model_name});
elseif length(authors) == 1
    lgd = legend(h,{model_name,authors{1}});
elseif length(authors) == 2
    lgd = legend(h,{model_name,authors{1},authors{2}});
end
set(lgd,'Location','best','FontSize',plot_opt.axes_size*0.6);

end