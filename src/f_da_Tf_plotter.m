function f_da_Tf_plotter(interp,outputs,source,model,tabgp,plot_opt,p)

% Unpack
iti = interp.iti;
itt = interp.itt;
ittf = interp.ittf;
range = interp.range;
time_ref = interp.time_ref;
alpha_ref = interp.alpha_ref; if source == "OSU", alpha_ref = alpha_ref*pi/180; end
cn_ref = interp.cn_ref;
cm_ref = interp.cm_ref;
cc_ref = interp.cc_ref;

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
[t,alpha,alpha_plunge,alpha_QS,q_QS,qR,R,alpha_cr,theta,theta_min,theta_max,S,P,T,alpha1_n,alpha1_m,alpha1_c,dalpha1_n,dalpha1_m,dalpha1_c,f_n,f_m,f_c,fprime_n,fprime_m,fprime_c,Tf_n,Tf_m,Tf_c,Ta_theta,alpha_C,c_n,c_nC,c_nI,c_nf,c_nv,c_m,c_mC,c_mI,c_mf,c_mv,dCP,c_c,c_l,c_d,alpha_lag,f2prime_n,f2prime_m,f2prime_c,RD,RD_theta,RD_tv0,f_diff_tv0,TvL_tv0,f_diff_tv0_2] = unpack_outputs(outputs,model);

% Define limits
if any(diff(plot_opt.xlim_vec))
    XLim = plot_opt.xlim_vec;
else
    XLim = plot_opt.xlim_alpha_plunge_vec;
    alpha = alpha_plunge;
end

%% Separation points - Normal coefficient
% Find interpolated values
c_nf_exp = cn_ref-c_nv(range)-c_nI(range);        % Estimated value for c_nf from the experiment (assuming the model's c_nv and c_nI are true)
K_f_exp = c_nf_exp./c_nC(range);                  % Kirchhoff factor
K_f_exp((K_f_exp>1)) = 1;                         % Bound upper limit
K_f_exp((K_f_exp<0)) = 0;                         % Bound lower limit
f2prime_n_exp = (2*real(K_f_exp.^(1/2))-1).^2;    % Estimated value for f2prime_n from the experiment
% Plot
tab1 = uitab('Parent',tabgp,'Title','f_n','BackgroundColor',[1 1 1]);
ax1 = axes('Parent',tab1,'FontSize',plot_opt.axes_size,'FontName','times new roman');
ax1.XLim = XLim; ax1.YLim = [0 1]; hold(ax1,'on'); 
plot(alpha(iti:end)*180/pi,f_n(iti:end),'k-',alpha(iti:end)*180/pi,fprime_n(iti:end),'b-',alpha(iti:end)*180/pi,f2prime_n(iti:end),'r-',alpha_ref*180/pi,f2prime_n_exp,'m-',alpha(iti:end)*180/pi,qR(iti:end),'c-',alpha(iti:end)*180/pi,RD(iti:end),'g-','LineWidth',plot_opt.lw,'Parent',ax1);
xlabel('$\alpha$ [deg]','FontWeight','normal','FontSize',plot_opt.axes_size);
ylabel('Delayed separation points','FontWeight','normal','FontSize',plot_opt.axes_size);
legend('f_n','f^{\prime}_n','f^{\prime\prime}_n','f^{\prime\prime}_n exp.','qR','R^{\prime}');
set(legend,'Location','best','FontSize',plot_opt.axes_size*0.6);
grid on

%% Separation points - Moment coefficient
% Find interpolated values
c_mf_interp = cm_ref-c_mv(range)-c_mI(range); % Estimated value for c_mf from the experiment (assuming the model's c_mv and c_mI are true)
% Newton-Raphson iterations to find f2prime_m
f2prime_m_interp = nan(ittf-itt+1,1);
if max(RD) > 0.05
    for i=1:ittf-itt+1
        if i == 1, f_guess = 0.5; else, f_guess = f2prime_m_interp(i-1); end
        j = i-1+itt;
        f2prime_m_interp(i) = fsolve(@(x) find_f2prime_m(x,c_mf_interp(i),c_nf(j),R(j),RD(j),S(j),theta(j),q_QS(j),p.kappa_0,p.kappa_1,p.kappa_2,p.kappa_3,p.K0,p.K1,p.K2),f_guess,optimoptions('fsolve','Display','off')); % Solve transcendental equation for f2prime_m
        f2prime_m_interp(i) = real(max([0,min([1,f2prime_m_interp(i)])]));
    end
end
% Plot
tab2 = uitab('Parent',tabgp,'Title','f_m','BackgroundColor',[1 1 1]);
ax2 = axes('Parent',tab2,'FontSize',plot_opt.axes_size,'FontName','times new roman');
ax2.XLim = XLim; ax2.YLim = [0 1]; hold(ax2,'on'); 
plot(alpha(iti:end)*180/pi,f_m(iti:end),'k-',alpha(iti:end)*180/pi,fprime_m(iti:end),'b-',alpha(iti:end)*180/pi,f2prime_m(iti:end),'r-',alpha_ref*180/pi,real(f2prime_m_interp),'m-',alpha(iti:end)*180/pi,qR(iti:end),'c-',alpha(iti:end)*180/pi,RD(iti:end),'g-','LineWidth',plot_opt.lw,'Parent',ax2);
xlabel('$\alpha$ [deg]','FontWeight','normal','FontSize',plot_opt.axes_size);
ylabel('Delayed separation points','FontWeight','normal','FontSize',plot_opt.axes_size);
legend('f_m','f^{\prime}_m','f^{\prime\prime}_m',' f^{\prime\prime}_m exp.','qR','R^{\prime}');
set(legend,'Location','best','FontSize',plot_opt.axes_size*0.6);
grid on

%% Separation points - Chordwise coefficient
f2prime_c_interp = nan(ittf-itt+1,1);
theta_c = abs(alpha_lag)/p.alpha1_0c;
for i=1:ittf-itt+1
    upstroke = theta(i-1+itt)*q_QS(i-1+itt) >= 0;
    if theta_c(i-1+itt) < 1
        f2prime_c_interp(i) = real(((cc_ref(i)+p.c_d0*cos(alpha_ref(i)))/(((1-p.eta*(RD(i-1+itt)^2+p.E0*S(i-1+itt)*R(i-1+itt)*theta(i-1+itt)^(1/2)*~upstroke)))*(p.c_n_alpha*sin((alpha_C(i-1+itt)-p.alpha_0L))^2)))^(1/(1/2+abs(theta_c(i-1+itt)))));
    else
        f2prime_c_interp(i) = real(((cc_ref(i)-p.E1*min([1,theta_c(i-1+itt)^3-1])+p.c_d0*cos(alpha_ref(i)))/(((1-p.eta*(RD(i-1+itt)^2+p.E0*S(i-1+itt)*R(i-1+itt)*theta(i-1+itt)^(1/2)*~upstroke)))*(p.c_n_alpha*sin((alpha_C(i-1+itt)-p.alpha_0L))^2)))^(1/(1/2+abs(theta_c(i-1+itt)))));
    end
end
tab3 = uitab('Parent',tabgp,'Title','f_c','BackgroundColor',[1 1 1]);
ax3 = axes('Parent',tab3,'FontSize',plot_opt.axes_size,'FontName','times new roman');
ax3.XLim = XLim; ax3.YLim = [0 1]; hold(ax3,'on'); 
% Plot
plot(alpha(iti:end)*180/pi,f_c(iti:end),'k-',alpha(iti:end)*180/pi,fprime_c(iti:end),'b-',alpha(iti:end)*180/pi,f2prime_c(iti:end),'r-',alpha_ref*180/pi,f2prime_c_interp,'m-',alpha(iti:end)*180/pi,qR(iti:end),'c-',alpha(iti:end)*180/pi,RD(iti:end),'g-','LineWidth',plot_opt.lw,'Parent',ax3);
xlabel('$\alpha$ [deg]','FontWeight','normal','FontSize',plot_opt.axes_size);
ylabel('Delayed separation points','FontWeight','normal','FontSize',plot_opt.axes_size);
legend('f_c','f^{\prime}_c','f^{\prime\prime}_c','f^{\prime\prime}_c exp.','qR','R^{\prime}');
set(legend,'Location','best','FontSize',plot_opt.axes_size*0.5);
grid on

%% Angle offsets
tab4 = uitab('Parent',tabgp,'Title','d_a1','BackgroundColor',[1 1 1]);
ax4 = axes('Parent',tab4,'FontSize',plot_opt.axes_size,'FontName','times new roman');
ax4.XLim = XLim;
hold(ax4,'on'); 
plot(alpha(iti:end)*180/pi,dalpha1_n(iti:end)*180/pi,'k-',alpha(iti:end)*180/pi,dalpha1_m(iti:end)*180/pi,'b-',alpha(iti:end)*180/pi,dalpha1_c(iti:end)*180/pi,'g-','LineWidth',plot_opt.lw,'Parent',ax4);
xlabel('$\alpha$ [deg]','FontWeight','normal','FontSize',plot_opt.axes_size);
ylabel('$\delta_{\alpha_1}$ [deg]','FontWeight','normal','FontSize',plot_opt.axes_size);
legend('Normal','Moment','Chordwise');
set(legend,'Location','best','FontSize',plot_opt.axes_size*0.6);
grid on

%% Time delay constants
tab5 = uitab('Parent',tabgp,'Title','T_f','BackgroundColor',[1 1 1]);
ax5 = axes('Parent',tab5,'FontSize',plot_opt.axes_size,'FontName','times new roman');
ax5.XLim = XLim;
hold(ax5,'on'); 
plot(alpha(iti:end)*180/pi,Tf_n(iti:end)/p.Tf,'k-',alpha(iti:end)*180/pi,Tf_m(iti:end)/p.Tf,'b-',alpha(iti:end)*180/pi,Tf_c(iti:end)/p.Tf,'g-','LineWidth',plot_opt.lw,'Parent',ax5);
xlabel('$\alpha$ [deg]','FontWeight','normal','FontSize',plot_opt.axes_size);
ylabel('$T_f/T_{f_0}$','FontWeight','normal','FontSize',plot_opt.axes_size);
legend('Normal','Moment','Chordwise');
set(legend,'Location','best','FontSize',plot_opt.axes_size*0.6);
grid on

%% Stall onset ratio
tab6 = uitab('Parent',tabgp,'Title','theta','BackgroundColor',[1 1 1]);
ax6 = axes('Parent',tab6,'FontSize',plot_opt.axes_size,'FontName','times new roman');
ax6.XLim = XLim;
hold(ax6,'on'); 
plot(alpha(iti:end)*180/pi,theta(iti:end),'k-',alpha(iti:end)*180/pi,R(iti:end),'b-',alpha(iti:end)*180/pi,RD(iti:end),'r-',alpha(iti:end)*180/pi,RD(iti:end).^p.lambda_2,'m-',alpha(iti:end)*180/pi,RD_theta(iti:end),'g-','LineWidth',plot_opt.lw,'Parent',ax6);
xlabel('$\alpha$ [deg]','FontWeight','normal','FontSize',plot_opt.axes_size);
ylabel('$\theta$','FontWeight','normal','FontSize',plot_opt.axes_size);
R_prime_to_lambda_2 = "R^{\prime}^{" + num2str(p.lambda_2) + "}";
legend('\theta','R','R^{\prime}',R_prime_to_lambda_2,'R^{\prime}_{\theta}');
set(legend,'Location','best','FontSize',plot_opt.axes_size*0.6);
grid on

end