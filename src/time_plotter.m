function [h,fig] = time_plotter(tabbed,fig,tab_name,coef_name,t,coef,time_interp,coef_interp,time_exp,coef_exp,time_mod,coef_mod,time_cycles,coef_cycles,model,outputs,params,interp,authors,source,tabgp,model_name,plot_opt)

% Set plot style for general harmonic cases
if contains(source,"hargen_test")
    switch params.time_plot_style
        case "pitch"
            source = "other";
        case "gust"
            source = "gust_test";
        case "tvfs"
            source = "tvfs_test";
    end
end

% Unpack and set arrays
switch source
    case {"NASA","GU","OSU","other","flap_test"}
        t_cycle = params.t_cycle; itt = interp.itt; ittf = interp.ittf;
        t_plot = (t(itt:ittf)-t(itt))/t_cycle*360-90;
        coef_plot = coef(itt:ittf);
    case "gust_test"
        U = params.U; b = params.b; alpha_g = atan(params.gust_options.wg_0/params.U); c_n_alpha = params.c_n_alpha;
        plot_increment = params.plot_increment; 
        gust_profile = params.gust_profile;
        % Set plot mode
        t_plot = U/b*t;
        if plot_increment
            coef_plot = coef-coef(1);
            coef_exp = coef_exp-coef_exp(1);
            coef_mod = coef_mod-coef_mod(1);
        else
            if gust_profile == "sharp-edge"
                if ismember(tab_name,{'cl x t','cn x t','cm x t'})
                    coef_plot = coef/c_n_alpha/alpha_g;
                else
                    coef_plot = coef;
                end
            else
                coef_plot = coef;
            end
        end
    case "tvfs_test"
        t_cycle = params.t_cycle; itt = interp.itt; ittf = interp.ittf;
        % Plot lift and normal coefficients' ratios to quasi-steady values
        t_plot = (t(itt:ittf)-t(itt))/t_cycle*360;
        coef_plot = coef(itt:ittf);
        if ismember(tab_name,{'cl x t','cn x t'})
            coefs_qs = params.c_n_alpha*outputs.alpha(1);
            coef_plot = coef_plot./coefs_qs;
        end
    case "hargen_test"
        itt = interp.itt; ittf = interp.ittf;
        t_plot = t(itt:ittf);
        coef_plot = coef(itt:ittf);
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
ax.XLabel.FontWeight = 'normal'; ax.XLabel.FontSize = plot_opt.axes_size;
ax.YLabel.FontWeight = 'normal'; ax.YLabel.FontSize = plot_opt.axes_size;
ax.Title.String = plot_opt.title_str; ax.Title.FontWeight = 'normal'; ax.Title.FontSize = plot_opt.axes_size*0.6;
switch source
    case "gust_test"
        ax.XLim = U/b*[t(1) t(end)];
        ax.XLabel.String = "$\tau$ [semi-chords]";
        if plot_increment
            ax.YLabel.String = coef_name + " increment";
        else
            if gust_profile == "sharp-edge"
                if ismember(tab_name,{'cl x t','cn x t','cm x t'})
                    ax.YLabel.String = "Normalized " + coef_name;
                else
                    ax.YLabel.String = coef_name;
                end
            else
                ax.YLabel.String = coef_name;
            end
        end
    case "tvfs_test"
        ax.XLim = [0 360]; ax.XTick = 0:45:360;
        ax.XLabel.String = '$\omega t$ [deg]';
        if ismember(tab_name,{'cl x t','cn x t'})
            coef_name = coef_name + " ratio";
        end
        ax.YLabel.String = coef_name;
    case "hargen_test"
        ax.XLabel.String = '$t$ [s]';
        ax.XLabel.String = coef_name;
    otherwise
        ax.XLim = [-90 270]; ax.XTick = -90:90:270;
        ax.XLabel.String = '$\omega t$ [deg]';
        ax.YLabel.String = coef_name;
end

%% Plot
switch source
    case "NASA"
        if ismember(tab_name,{'cl x t','cm x t','cd x t'})
            h = plot(t_plot,coef_plot,'k-',time_exp,coef_exp,'ko','LineWidth',plot_opt.lw,'MarkerSize',plot_opt.ms,'Parent',ax);
        else
            h = plot(t_plot,coef_plot,'k-',time_interp,coef_interp,'b-','LineWidth',plot_opt.lw,'Parent',ax);
        end
    case "GU"
        h = plot(t_plot,coef_plot,'k-',time_exp,coef_exp,'ko','LineWidth',plot_opt.lw,'MarkerSize',plot_opt.ms,'Parent',ax);
    case "OSU"
        % Plot model results
        h_mod = plot(t_plot,coef_plot,'k-','LineWidth',plot_opt.lw,'Parent',ax);
        % Plot the multiple experimental cycles
        for n=1:length(time_cycles)
            plot(time_cycles{n}(1:end-1),coef_cycles{n}(1:end-1),'bo','LineWidth',plot_opt.lw,'MarkerSize',plot_opt.ms,'Parent',ax);
        end
        % Plot mean interpolated values
        h_exp = plot(time_interp,coef_interp,'b-','LineWidth',plot_opt.lw,'MarkerSize',plot_opt.ms,'Parent',ax);
        % Legend
        h = plot(nan,nan,'k-',nan,nan,'b-o','LineWidth',plot_opt.lw,'MarkerSize',plot_opt.ms,'Parent',ax);
    otherwise
        if isnan(time_exp(1)) && isempty(authors{1})
            h = plot(t_plot,coef_plot,'k-','LineWidth',plot_opt.lw,'MarkerSize',plot_opt.ms,'Parent',ax);
        elseif length(authors) == 1  || isnan(time_mod(1)) 
            h = plot(t_plot,coef_plot,'k-',time_exp,coef_exp,'ko','LineWidth',plot_opt.lw,'MarkerSize',plot_opt.ms,'Parent',ax);
        elseif length(authors) == 2
            h = plot(t_plot,coef_plot,'k-',time_exp,coef_exp,'ko',time_mod,coef_mod,'k--','LineWidth',plot_opt.lw,'MarkerSize',plot_opt.ms,'Parent',ax);
        end            
end

% Set legend
if isnan(time_exp(1)) && isempty(authors{1})
    lgd = legend(h,{model_name});
elseif length(authors) == 1 || isnan(time_mod(1)) 
    lgd = legend(h,{model_name,authors{1}});
elseif length(authors) == 2
    lgd = legend(h,{model_name,authors{1},authors{2}});
end
set(lgd,'Location','best','FontSize',plot_opt.axes_size*0.6);

% Set correct handle output for OSU case
if source == "OSU", h = [h_mod; h_exp]; end

end