function [h,fig] = alpha_plotter(tabbed,fig,tab_name,coef_name,alpha,coef,alpha_interp,coef_interp,alpha_exp,coef_exp,alpha_mod,coef_mod,alpha_cycles,coef_cycles,authors,iti,source,tabgp,model_name,plot_opt)

% Check input
if diff(plot_opt.xlim_vec) == 0
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
ax.XLim = plot_opt.xlim_vec;
ax.XLabel.String = '$\alpha$ [deg]'; ax.XLabel.FontWeight = 'normal'; ax.XLabel.FontSize = plot_opt.axes_size;
ax.YLabel.String = coef_name; ax.YLabel.FontWeight = 'normal'; ax.YLabel.FontSize = plot_opt.axes_size;

%% Plot
alpha_plot = alpha(iti:end)*180/pi; coef_plot = coef(iti:end);
switch source
    case "NASA"
        if ismember(tab_name,{'cl','cm','cd'})
            h = plot(alpha_plot,coef_plot,'k-',alpha_exp,coef_exp,'ko','LineWidth',plot_opt.lw,'Parent',ax);
        else
            h = plot(alpha_plot,coef_plot,'k-',alpha_interp*180/pi,coef_interp,'b-','LineWidth',plot_opt.lw,'Parent',ax);
        end
    case "GU"
        h = plot(alpha_plot,coef_plot,'k-',alpha_exp,coef_exp,'ko','LineWidth',plot_opt.lw,'Parent',ax);
    case "OSU"
        % Plot model results
        h_mod = plot(alpha_plot,coef_plot,'k-','LineWidth',plot_opt.lw,'Parent',ax);
        % Plot the multiple experimental cycles
        for n=1:length(alpha_cycles)
            plot(alpha_cycles{n},coef_cycles{n},'bo','LineWidth',plot_opt.lw,'MarkerSize',plot_opt.ms,'Parent',ax);
        end
        % Plot mean interpolated values
        h_exp = plot(alpha_interp,coef_interp,'b-','LineWidth',plot_opt.lw,'MarkerSize',plot_opt.ms,'Parent',ax);
        % For legend
        h = plot(nan,nan,'k-',nan,nan,'b-o','LineWidth',plot_opt.lw,'MarkerSize',plot_opt.ms,'Parent',ax);
    otherwise
        if isnan(alpha_exp(1)) && isempty(authors{1})
            h = plot(alpha_plot,coef_plot,'k-','LineWidth',plot_opt.lw,'MarkerSize',plot_opt.ms,'Parent',ax);
        elseif length(authors) == 1 || isnan(alpha_mod(1))   
            h = plot(alpha_plot,coef_plot,'k-',alpha_exp,coef_exp,'ko','LineWidth',plot_opt.lw,'MarkerSize',plot_opt.ms,'Parent',ax);
        elseif length(authors) == 2
            h = plot(alpha_plot,coef_plot,'k-',alpha_exp,coef_exp,'ko',alpha_mod,coef_mod,'k--','LineWidth',plot_opt.lw,'MarkerSize',plot_opt.ms,'Parent',ax);
        end
end

% Set legend
if isnan(alpha_exp(1)) && isempty(authors{1})
    lgd = legend(h,{model_name});
elseif length(authors) == 1 || isnan(alpha_mod(1))
    lgd = legend(h,{model_name,authors{1}});
elseif length(authors) == 2
    lgd = legend(h,{model_name,authors{1},authors{2}});
end
set(lgd,'Location','best','FontSize',plot_opt.axes_size*0.6);

% Set correct handle output for OSU case
if source == "OSU", h = [h_mod; h_exp]; end

end