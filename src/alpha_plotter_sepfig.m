function alpha_plotter_sepfig(fig,coef_name,alpha,coef,alpha_interp,coef_interp,alpha_exp,coef_exp,alpha_mod,coef_mod,alpha_cycles,coef_cycles,data,params,interp,source,model_name,plot_opt,case_now)

% Unpack 
authors = data.authors;
M  = params.M;
k = params.k;
a_0 = params.a_0;
a_1 = params.a_1;
iti = interp.iti;

%% Initialize tab and set axes design
ax = axes('Parent',fig,'FontSize',plot_opt.axes_size,'FontName','times new roman');
ax.XLim = plot_opt.xlim_vec;
hold(ax,'on'); 
xlabel('\alpha [deg]','FontWeight','normal','FontSize',plot_opt.axes_size);
ylabel(coef_name,'FontWeight','normal','FontSize',plot_opt.axes_size);
grid on

%% Configure title
info_str = "M = " +  num2str(M,'%.3f') + ", k = " + num2str(k,'%.3f') + ", A0 = " +  num2str(a_0*180/pi,'%.1f') + "^{\circ}, A1 = " + num2str(a_1*180/pi,'%.1f') + "^{\circ}";
if length(authors) == 1
    if source == "NASA" 
        source_str = "Frame " + num2str(case_now) + ": ";
    elseif source == "GU"
        source_str = "GUD " + num2str(case_now) + ": ";
    elseif source == "OSU" 
        source_str = "OSU run " + num2str(case_now) + ": ";
    else
        source_str = "";
    end
else
    source_str = "";
end
title_str = source_str + info_str;
title(title_str,'FontWeight','normal','FontSize',plot_opt.axes_size*0.6);

%% Plot
if length(authors) == 1 
    if source == "NASA" 
        if ismember(coef_name,{'Lift coefficient','Moment coefficient','Chordwise coefficient'})
            plot(alpha(iti:end)*180/pi,coef(iti:end),'k-',alpha_exp,coef_exp,'ko','LineWidth',plot_opt.lw,'Parent',ax);
        else
            plot(alpha(iti:end)*180/pi,coef(iti:end),'k-',alpha_interp*180/pi,coef_interp,'b-','LineWidth',plot_opt.lw,'Parent',ax);
        end
        legend(model_name,authors{1});
    elseif source == "GU"
        plot(alpha(iti:end)*180/pi,coef(iti:end),'k-',alpha_exp,coef_exp,'ko','LineWidth',plot_opt.lw,'Parent',ax);
        legend(model_name,authors{1});    
    elseif source == "OSU"
        % Plot model results
        plot(alpha(iti:end)*180/pi,coef(iti:end),'k-','LineWidth',plot_opt.lw,'Parent',ax);
        % Plot the multiple experimental cycles
        n_cycles = length(alpha_cycles); 
        for n=1:n_cycles
            plot(alpha_cycles{n},coef_cycles{n},'bo','LineWidth',plot_opt.lw,'MarkerSize',plot_opt.ms,'Parent',ax);
        end
        % Plot mean interpolated values
        plot(alpha_interp,coef_interp,'b-','LineWidth',plot_opt.lw,'MarkerSize',plot_opt.ms,'Parent',ax);
        % Legend
        dummy_plot = plot(nan,nan,'k-',nan,nan,'b-o','LineWidth',plot_opt.lw,'MarkerSize',plot_opt.ms,'Parent',ax);
        legend(dummy_plot,{model_name,authors{1}});
    elseif source == "test"
        plot(alpha(iti:end)*180/pi,coef(iti:end),'k-','LineWidth',plot_opt.lw,'Parent',ax);
        legend(model_name);
    else
        plot(alpha(iti:end)*180/pi,coef(iti:end),'k-',alpha_exp,coef_exp,'ko','LineWidth',plot_opt.lw,'MarkerSize',plot_opt.ms,'Parent',ax);
        legend(model_name,authors{1});
    end
else
    plot(alpha(iti:end)*180/pi,coef(iti:end),'k-',alpha_exp,coef_exp,'ko',alpha_mod,coef_mod,'k--','LineWidth',plot_opt.lw,'MarkerSize',plot_opt.ms,'Parent',ax);
    if isnan(alpha_exp(1))
        legend(model_name);
    else
        legend(model_name,authors{1},authors{2});
    end
end

% Set legend
set(legend,'Location','best','FontSize',plot_opt.axes_size*0.6);

end