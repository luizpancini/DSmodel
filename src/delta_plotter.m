function h = delta_plotter(coef_name,tab_name,delta,coef,delta_exp,coef_exp,delta_mod,coef_mod,authors,iti,tabgp,model_name,plot_opt)

% Check input
if diff(plot_opt.xlim_delta_vec) == 0
    h = gobjects(0); return;
end

%% Initialize tab and set axes design
tab = uitab('Parent',tabgp,'Title',tab_name,'BackgroundColor',[1 1 1]);
ax = axes('Parent',tab,'FontSize',plot_opt.axes_size,'FontName','times new roman');
ax.XLim = plot_opt.xlim_delta_vec;
hold(ax,'on'); 
xlabel('$\delta$ [deg]','FontWeight','normal','FontSize',plot_opt.axes_size);
ylabel(coef_name,'FontWeight','normal','FontSize',plot_opt.axes_size);
title(plot_opt.title_str,'FontWeight','normal','FontSize',plot_opt.axes_size*0.6);
grid on

%% Plot
delta_plot = delta(iti:end)*180/pi; coef_plot = coef(iti:end);
if isnan(delta_exp)
    h = plot(delta_plot,coef_plot,'k-','LineWidth',plot_opt.lw,'MarkerSize',plot_opt.ms,'Parent',ax);
elseif length(authors) == 1
    h = plot(delta_plot,coef_plot,'k-',delta_exp,coef_exp,'ko','LineWidth',plot_opt.lw,'MarkerSize',plot_opt.ms,'Parent',ax);
elseif length(authors) == 2
    h = plot(delta_plot,coef_plot,'k-',delta_exp,coef_exp,'ko',delta_mod,coef_mod,'k--','LineWidth',plot_opt.lw,'MarkerSize',plot_opt.ms,'Parent',ax);
end
% Set legend
if isnan(delta_exp)
    legend(h,{model_name});
elseif length(authors) == 1
    legend(h,{model_name,authors{1}});
elseif length(authors) == 2
    legend(h,{model_name,authors{1},authors{2}});
end
set(legend,'Location','best','FontSize',plot_opt.axes_size*0.6);

end