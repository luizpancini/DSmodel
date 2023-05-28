function [h,plot_opt] = init_plotter(h,base_name,case_now,source,params,tabbed)

% Unpack
a_0 = params.a_0;
a_1 = params.a_1;
M = params.M;
k = params.k;
if isnumeric(case_now), case_str =  num2str(case_now); else, case_str = ""; end

% Default plot options
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultLegendInterpreter','tex')
axes_size = 20;
lw = 1;
ms = 5;

% Create figure and tab group
h.fig = figure('InvertHardcopy','off','Color',[1 1 1]); 
h.fig.Name = base_name + case_str;
h.fig.NumberTitle = 'off';
h.fig.PaperSize = [9 7.2];
h.fig.Position = [3.5540e+02 129 920 6.2960e+02];
if tabbed, h.tabgp = uitabgroup(h.fig); end

% Set axes limits
xlim_vec = sort(180/pi*[a_0-a_1, a_0+a_1]);
if isfield(params,'d_0')
    d_0 = params.d_0; d_1 = params.d_1;
    xlim_delta_vec = sort(180/pi*[d_0-d_1, d_0+d_1]);
else
    xlim_delta_vec = [];
end
if isfield(params,'a_1h')
    a_0 = params.a_0; a_1h = params.a_1h; k_h = params.k_h;
    xlim_alpha_plunge_vec = sort(180/pi*[a_0-a_1h, a_0+a_1h]);
else
    xlim_alpha_plunge_vec = [];
end

% Configure title
info_str = "";
switch source
    case {'NASA','GU','OSU','other'}
        info_str = ": $M$ = " +  num2str(M,'%.3f') + ", $k$ = " + num2str(k,'%.3f') + ", $\alpha_0$ = " +  num2str(a_0*180/pi,'%.1f') + "$^{\circ}$, $\alpha_1$ = " + num2str(a_1*180/pi,'%.1f') + "$^{\circ}$, $k_h$ = " + num2str(k_h,'%.3f') + ", $\alpha_{1_h}$ = " + num2str(a_1h*180/pi,'%.1f') + "$^{\circ}$";
    case 'gust_test'
        gust_profile = params.gust_profile; wg_0 = params.gust_options.wg_0; H = params.gust_options.H; U = params.U; lambda_g = params.lambda_g;
        alpha_g = wg_0/U*180/pi; 
        info_str = ": Gust profile: " + gust_profile + ", ";
        switch gust_profile
            case {"sharp-edge","1-cos"}
                info_str = info_str + "$\alpha_g$ = " + num2str(alpha_g,'%.1f') + "$^{\circ}$, $H$ = " + num2str(H,'%.0f');
            case "sine"
                f_sin = params.gust_options.f_sin;
                info_str = info_str + "$\alpha_g$ = " + num2str(alpha_g,'%.1f') + "$^{\circ}$, $H$ = " + num2str(H,'%.0f') + ", $f$ = " + num2str(f_sin,'%.2f') + " Hz";
            case "swept-sine"
                f_i = params.gust_options.f_i; f_f = params.gust_options.f_f;
                info_str = info_str + "$\alpha_g$ = " + num2str(alpha_g,'%.1f') + "$^{\circ}$, $H$ = " + num2str(H,'%.0f') + ", $f_i$ = " + num2str(f_i,'%.2f') + " Hz, $f_f$ = " + num2str(f_f,'%.2f') + " Hz";
        end          
        info_str = info_str + ", $M$ = " +  num2str(M,'%.3f') + ", $\alpha_0$ = " +  num2str(a_0*180/pi,'%.1f') + "$^{\circ}$" + ", $\lambda_g$ = " + num2str(lambda_g,'%.2f');
    case 'flap_test'
        d_0 = params.d_0; d_1 = params.d_1; k_f = params.k_f;
        info_str = ": $M$ = " +  num2str(M,'%.3f') + ", $k$ = " + num2str(k,'%.2f') + ", $\alpha_0$ = " +  num2str(a_0*180/pi,'%.1f') + "$^{\circ}$, $\alpha_1$ = " + num2str(a_1*180/pi,'%.1f') + "$^{\circ}$" + ", $k_f$ = " + num2str(k_f,'%.2f') + ", $\delta_0$ = " +  num2str(d_0*180/pi,'%.1f') + "$^{\circ}$, $\delta_1$ = " + num2str(d_1*180/pi,'%.1f') + "$^{\circ}$";
    case 'tvfs_test'
        U_0 = params.U_0; U_1 = params.U_1; k_U = params.k_U; if isfield(params,'psi_a'), psi = params.psi_a; else, psi = params.psi_Ua; end
        info_str = ": $M_0$ = " +  num2str(M,'%.3f') + ", $\lambda_U$ = " + num2str(U_1/U_0,'%.2f') + ", $k_U$ = " + num2str(k_U,'%.2f') + ", $\alpha_0$ = " +  num2str(a_0*180/pi,'%.1f') + "$^{\circ}$, $\alpha_1$ = " + num2str(a_1*180/pi,'%.1f') + "$^{\circ}$" + ", $k$ = " + num2str(k,'%.2f') + ", $\psi$ = " + num2str(rad2deg(psi),'%.0f') + "$^{\circ}$";
end
title_str = base_name + case_str + info_str;

% Save to struct
plot_opt = variables2struct(struct,axes_size,lw,ms,xlim_vec,xlim_delta_vec,xlim_alpha_plunge_vec,title_str);

end