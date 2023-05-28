
% Initialize figures
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultLegendInterpreter','tex')
axes_size = 20;
lw = 1;
ms = 5;

% Create figure, axes and tab group
figure1 = figure('InvertHardcopy','off','Color',[1 1 1]); 
figure1.Position = [3.5540e+02 129 920 6.2960e+02];
tabgp = uitabgroup(figure1);
tab1 = uitab('Parent',tabgp,'Title','c_n','BackgroundColor',[1 1 1]); tab2 = uitab('Parent',tabgp,'Title','c_m','BackgroundColor',[1 1 1]);
ax1 = axes('Parent',tab1,'FontSize',axes_size,'FontName','times new roman');
ax1.XLim = [-90,270]; ax1.XLabel.String = "Cycle angle [deg]"; ax1.YLabel.String = "$c_n^v$"; hold(ax1,'on');
ax2 = axes('Parent',tab2,'FontSize',axes_size,'FontName','times new roman');
ax2.XLim = [-90,270]; ax2.XLabel.String = "Cycle angle [deg]"; ax2.YLabel.String = "$c_m^v$"; hold(ax2,'on');

% Colors 
c = lines(length(OUTPUTS));

% Loop over cases
for i=1:length(OUTPUTS)
    
    % Unpack interpolation data
    iti = OUTPUTS(i).interp.iti;
    itt = OUTPUTS(i).interp.itt;
    ittf = OUTPUTS(i).interp.ittf;
    range = OUTPUTS(i).interp.range;
    time_ref = OUTPUTS(i).interp.time_ref;
    alpha_ref = OUTPUTS(i).interp.alpha_ref; if INPUTS.source == "OSU", alpha_ref = alpha_ref*pi/180; end
    cn_ref = OUTPUTS(i).interp.cn_ref;
    cm_ref = OUTPUTS(i).interp.cm_ref;
    cc_ref = OUTPUTS(i).interp.cc_ref;
    
    % Re-interpolate reference values in the range of the last cycle
    N_range = length(range); interp_method = 'pchip';
    time_interp = linspace(-90,270,N_range);
    alpha_ref = interp1(time_ref,alpha_ref,time_interp,interp_method);
    cn_ref = interp1(time_ref,cn_ref,time_interp,interp_method);
    cm_ref = interp1(time_ref,cm_ref,time_interp,interp_method);
    if ~isnan(cc_ref(1))
        cc_ref = interp1(time_ref,cc_ref,time_interp,interp_method);
    else
        cc_ref = nan*cn_ref;
    end
    
    % Unpack output variables
    [t,alpha,c_n,c_m,c_c,c_l,c_d,c_nf,c_nI,c_mf,c_mI,c_nv,c_mv,c_cv,c_lv,f,tau_v,dCP,so_lim,qR,alpha1_n,K_f,c_nC,Tf_n,dalpha1_n,fprime_n,fprime_m,fprime_c,q,dalpha1_m,dalpha1_c,R_dot,Tf_m,Tf_c,theta,theta_min,theta_max,P,S,alpha_E,R,so_state,f2prime_n,f2prime_m,f2prime_c,RD,RD_theta,alpha_plunge] = unpack_outputs(OUTPUTS(i).outputs,INPUTS.model);
    t_range = (t(range)-t(range(1)))/OUTPUTS(i).params.t_cycle*360-90; alpha = alpha(range);
    c_n = c_n(range); c_nf = c_nf(range); c_nI = c_nI(range); c_nv = c_nv(range);
    c_m = c_m(range); c_mf = c_mf(range); c_mI = c_mI(range); c_mv = c_mv(range);
    
    % Estimated experimental data's vortex loads
    c_nv_est = cn_ref-(c_nf+c_nI);
    c_mv_est = cm_ref-(c_mf+c_mI);
    
    % Compare coefficients
    plot(t_range,c_nv,'-',time_interp,c_nv_est,'--','Color',c(i,:),'LineWidth',lw,'Parent',ax1);
    plot(t_range,c_mv,'-',time_interp,c_mv_est,'--','Color',c(i,:),'LineWidth',lw,'Parent',ax2);
    
end
lgd = legend('Model','Exp. (estimated)');
set(lgd,'Location','best','FontSize',12);