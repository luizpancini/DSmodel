% Compares model with my modifications and that of Sheng's.
clc
clear 
close all

%% Time span
n_discard = 3;              % Number of cycles to discard for plots
n_cycles = n_discard+2;     % Number of cycles to run
    
%% Choose experimental data to load from McAlister's frame, GU's data, some other experimental run or a desired test condition
GUD = 11012812;
frame = 10212;
run = 10;
test = {0.072; 0.248; 9.9*pi/180; 10.7*pi/180; 'NACA0012'}; % [M; k; a_0; a_1; airfoil]
INPUT = GUD;
[authors,data,U,beta,M,k,a_0,a_1,b,a_inf,airfoil,params] = read_data(INPUT);
lin_reat = 1; % Choose to use linear reattachment modelling for BLS

%% Indicial parameters
params = read_indicial_params(M,beta,b,a_inf,airfoil,params);

%% State space matrices
[A,B] = get_SS_matrices(U,b,beta,M,a_inf,airfoil,params);

%% Initial conditions: {x0} = 0
N = 14; % Number of states
x0 = zeros(N,1);
y0 = [a_0; 0; 0; 0];
y0(5:35) = nan;

%% ODE solver
dt = 1e-4;              % Maximum time step
t_cycle = 2*pi*b/(U*k); % Time of one cycle
tf = t_cycle*n_cycles;  % Test time
tspan = [0 tf];
options = {'hlim',dt};

%% BL modified model
tic; [t,x,y,xdot] = BL_RKF45(tspan,x0,y0,A,B,U,M,b,a_0,a_1,k,beta,params,options);
time_BL = toc;
% Output variables
[alpha,c_n,c_m,c_c,c_l,c_d,c_nf,c_nI,c_mf,c_mI,c_nv,...
c_mv,c_cv,c_lv,f,f2prime,tau_v,dCP,so_lim,qR,alpha1,...
K_f,c_nC,Tf_n,dalpha1,fprime,fprime_cm,fprime_cc,q,dalpha1_cm,...
dalpha1_cc,R_dot,Tf_m,Tf_c,theta,theta_min,theta_max,P,S,R,...
RD_theta,alpha_E] = BL_output_vars(x,y);
% Interpolate data and find normalized error coefficients
[iti,itt,ittf,time_interp,alpha_interp,cn_interp,cm_interp,cc_interp,cn_NRMSE,cm_NRMSE,cc_NRMSE] = interp_data(c_n,c_m,c_c,t,t_cycle,n_discard,data,authors,a_0,a_1,INPUT,0);

%% BLS model 
tic; [alpha_BLS,cn_BLS,cm_BLS,cc_BLS,cl_BLS,cd_BLS,t_BLS] = run_BLS(INPUT,lin_reat,dt,n_cycles);
time_BLS = toc;
% Interpolate data and find normalized error coefficients
[iti_BLS,itt_BLS,ittf_BLS,time_interp_BLS,alpha_interp_BLS,cn_interp_BLS,cm_interp_BLS,cc_interp_BLS,cn_BLS_NRMSE,cm_BLS_NRMSE,cc_BLS_NRMSE] = interp_data(cn_BLS,cm_BLS,cc_BLS,t_BLS,t_cycle,n_discard,data,authors,a_0,a_1,INPUT,0);

%% Compare models
cne_models = [cn_NRMSE; cn_BLS_NRMSE];
cme_models = [cm_NRMSE; cm_BLS_NRMSE];
cce_models = [cc_NRMSE; cc_BLS_NRMSE];
runtime_models = [time_BL; time_BLS];
comparison_table = table(cne_models,cme_models,cce_models,runtime_models); 
comparison_table.Properties.RowNames = {'BL mod','BLS'}; comparison_table.Properties.VariableNames = {'cn_error','cm_error','cc_error','runtime'};
comparison_table
% Computational time comparison
disp(['Time ratio BL/BLS: ' num2str(time_BL/time_BLS,'%.2f')]);

%% Plots
set(0,'DefaultTextInterpreter','tex')
set(0,'DefaultLegendInterpreter','tex')
axes_size = 20;
lw = 1;
ms = 5;
% Determine axes limits
if a_0+a_1 > a_0-a_1
    xlim_vec = [(a_0-a_1)*180/pi (a_0+a_1)*180/pi];
else
    xlim_vec = [(a_0+a_1)*180/pi (a_0-a_1)*180/pi];
end

% c_l
figure1 = figure('InvertHardcopy','off','Color',[1 1 1]); figure1.PaperSize = [6 4.5];
axes1 = axes('Parent',figure1,'FontSize',axes_size,'FontName','times new roman');
axes1.XLim = xlim_vec;
hold(axes1,'on'); 
if length(authors) == 1 
    if isnan(data{1}(1)) || isnan(authors{1}(1))
        plot1 = plot(alpha(iti:end)*180/pi,c_l(iti:end),'k-',alpha_BLS(iti_BLS:end)*180/pi,cl_BLS(iti_BLS:end),'k--','LineWidth',lw,'Parent',axes1);
        legend('Present model','Sheng et al. model');
    else
        plot1 = plot(alpha(iti:end)*180/pi,c_l(iti:end),'k-',alpha_BLS(iti_BLS:end)*180/pi,cl_BLS(iti_BLS:end),'k--',data{1,1},data{2,1},'ko','LineWidth',lw,'MarkerSize',ms,'Parent',axes1);
        legend('Present model','Sheng et al. model',authors{1});
    end
else
    plot1 = plot(alpha(iti:end)*180/pi,c_l(iti:end),'k-',data{1,1},data{2,1},'ko',data{17,1},data{18,1},'k--','LineWidth',lw,'MarkerSize',ms,'Parent',axes1);
    if isnan(data{1}(1))
        legend('Present model');
    else
        legend('Present model',authors{1},authors{2});
    end
end
xlabel('\alpha [deg]','FontWeight','normal','FontSize',axes_size);
ylabel('Lift coefficient','FontWeight','normal','FontSize',axes_size);
if strcmp(authors{1},'McAlister et al. (1982)') && length(authors) == 1
    title(['Frame ' num2str(frame) ': M = ' num2str(M) ', k = ' num2str(k) ', A0 = ' num2str(a_0*180/pi,'%.1f') '^{\circ}, A1 = ' num2str(a_1*180/pi,'%.1f') '^{\circ}'],'FontWeight','normal','FontSize',axes_size*0.6);
elseif strcmp(authors{1},'GU') && length(authors) == 1
    title(['GUD ' num2str(GUD) ': M = ' num2str(M,'%.3f') ', k = ' num2str(k,'%.2f') ', A0 = ' num2str(a_0*180/pi,'%.1f') '^{\circ}, A1 = ' num2str(a_1*180/pi,'%.1f') '^{\circ}'],'FontWeight','normal','FontSize',axes_size*0.6);
end
set(legend,'Location','best','FontSize',axes_size*0.6);
grid on

% c_m
figure2 = figure('InvertHardcopy','off','Color',[1 1 1]); figure2.PaperSize = [6 4.5];
axes2 = axes('Parent',figure2,'FontSize',axes_size,'FontName','times new roman');
axes2.XLim = xlim_vec;
hold(axes2,'on'); 
if length(authors) == 1 
    if isnan(data{3}(1)) || isnan(authors{1}(1))
        plot2 = plot(alpha(iti:end)*180/pi,c_m(iti:end),'k-',alpha_BLS(iti_BLS:end)*180/pi,cm_BLS(iti_BLS:end),'k--','LineWidth',lw,'Parent',axes2);
        legend('Present model','Sheng et al. model');
    else
        plot2 = plot(alpha(iti:end)*180/pi,c_m(iti:end),'k-',alpha_BLS(iti_BLS:end)*180/pi,cm_BLS(iti_BLS:end),'k--',data{3,1},data{4,1},'ko','LineWidth',lw,'MarkerSize',ms,'Parent',axes2);
        legend('Present model','Sheng et al. model',authors{1});
    end
else
    plot2 = plot(alpha(iti:end)*180/pi,c_m(iti:end),'k-',data{3,1},data{4,1},'ko',data{19,1},data{20,1},'k--','LineWidth',lw,'MarkerSize',ms,'Parent',axes2);
    if isnan(data{3}(1))
        legend('Present model');
    else
        legend('Present model',authors{1},authors{2});
    end
end
xlabel('\alpha [deg]','FontWeight','normal','FontSize',axes_size);
ylabel('Moment coefficient','FontWeight','normal','FontSize',axes_size);
if strcmp(authors{1},'McAlister et al. (1982)') && length(authors) == 1
    title(['Frame ' num2str(frame) ': M = ' num2str(M) ', k = ' num2str(k) ', A0 = ' num2str(a_0*180/pi,'%.1f') '^{\circ}, A1 = ' num2str(a_1*180/pi,'%.1f') '^{\circ}'],'FontWeight','normal','FontSize',axes_size*0.6);
elseif strcmp(authors{1},'GU') && length(authors) == 1
    title(['GUD ' num2str(GUD) ': M = ' num2str(M,'%.3f') ', k = ' num2str(k,'%.2f') ', A0 = ' num2str(a_0*180/pi,'%.1f') '^{\circ}, A1 = ' num2str(a_1*180/pi,'%.1f') '^{\circ}'],'FontWeight','normal','FontSize',axes_size*0.6);
end
set(legend,'Location','best','FontSize',axes_size*0.6);
grid on

% c_d
figure3 = figure('InvertHardcopy','off','Color',[1 1 1]); figure3.PaperSize = [6 4.5];
axes3 = axes('Parent',figure3,'FontSize',axes_size,'FontName','times new roman');
axes3.XLim = xlim_vec;
hold(axes3,'on'); 
if length(authors) == 1 
    if isnan(data{5}(1)) || isnan(authors{1}(1))
        plot3 = plot(alpha(iti:end)*180/pi,c_d(iti:end),'k-',alpha_BLS(iti_BLS:end)*180/pi,cd_BLS(iti_BLS:end),'k--','LineWidth',lw,'Parent',axes3);
        legend('Present model','Sheng et al. model');
    else
        plot3 = plot(alpha(iti:end)*180/pi,c_d(iti:end),'k-',alpha_BLS(iti_BLS:end)*180/pi,cd_BLS(iti_BLS:end),'k--',data{5,1},data{6,1},'ko','LineWidth',lw,'MarkerSize',ms,'Parent',axes3);
        legend('Present model','Sheng et al. model',authors{1});
    end
else
    plot3 = plot(alpha(iti:end)*180/pi,c_d(iti:end),'k-',data{5,1},data{6,1},'ko',data{21,1},data{22,1},'k--','LineWidth',lw,'MarkerSize',ms,'Parent',axes3);
    if isnan(data{5}(1))
        legend('Present model');
    else
        legend('Present model',authors{1},authors{2});
    end
end
xlabel('\alpha [deg]','FontWeight','normal','FontSize',axes_size);
ylabel('Drag coefficient','FontWeight','normal','FontSize',axes_size);
if strcmp(authors{1},'McAlister et al. (1982)') && length(authors) == 1
    title(['Frame ' num2str(frame) ': M = ' num2str(M) ', k = ' num2str(k) ', A0 = ' num2str(a_0*180/pi,'%.1f') '^{\circ}, A1 = ' num2str(a_1*180/pi,'%.1f') '^{\circ}'],'FontWeight','normal','FontSize',axes_size*0.6);
elseif strcmp(authors{1},'GU') && length(authors) == 1
    title(['GUD ' num2str(GUD) ': M = ' num2str(M,'%.3f') ', k = ' num2str(k,'%.2f') ', A0 = ' num2str(a_0*180/pi,'%.1f') '^{\circ}, A1 = ' num2str(a_1*180/pi,'%.1f') '^{\circ}'],'FontWeight','normal','FontSize',axes_size*0.6);
end
set(legend,'Location','best','FontSize',axes_size*0.6);
grid on

% c_n
figure8 = figure('InvertHardcopy','off','Color',[1 1 1]); figure8.PaperSize = [6 4.5];
axes8 = axes('Parent',figure8,'FontSize',axes_size,'FontName','times new roman');
axes8.XLim = xlim_vec;
hold(axes8,'on'); 
if length(authors) == 1 
    if ~isnan(alpha_interp(1)) && strcmp(authors{1},'McAlister et al. (1982)')
        plot8 = plot(alpha(iti:end)*180/pi,c_n(iti:end),'k-',alpha_BLS(iti_BLS:end)*180/pi,cn_BLS(iti_BLS:end),'k--',alpha_interp*180/pi,cn_interp,'b-','LineWidth',lw,'Parent',axes8);
        legend('Present model','Sheng et al. model',[authors{1}]);  
    elseif isnan(data{13}(1)) || isnan(authors{1}(1))
        plot8 = plot(alpha(iti:end)*180/pi,c_n(iti:end),'k-',alpha_BLS(iti_BLS:end)*180/pi,cn_BLS(iti_BLS:end),'k--','LineWidth',lw,'Parent',axes8);
        legend('Present model','Sheng et al. model');
    else
        plot8 = plot(alpha(iti:end)*180/pi,c_n(iti:end),'k-',alpha_BLS(iti_BLS:end)*180/pi,cn_BLS(iti_BLS:end),'k--',data{13,1},data{14,1},'ko','LineWidth',lw,'MarkerSize',ms,'Parent',axes8);
        legend('Present model','Sheng et al. model',authors{1});
    end
else
    plot8 = plot(alpha(iti:end)*180/pi,c_n(iti:end),'k-',data{13,1},data{14,1},'ko',data{23,1},data{24,1},'k--','LineWidth',lw,'MarkerSize',ms,'Parent',axes8);
    if isnan(data{13}(1))
        legend('Present model');
    else
        legend('Present model',authors{1},authors{2});
    end
end
set(legend,'Location','best','FontSize',axes_size*0.6);
xlabel('\alpha [deg]','FontWeight','normal','FontSize',axes_size);
ylabel('Normal coefficient','FontWeight','normal','FontSize',axes_size);
if strcmp(authors{1},'McAlister et al. (1982)') && length(authors) == 1
    title(['Frame ' num2str(frame) ': M = ' num2str(M) ', k = ' num2str(k) ', A0 = ' num2str(a_0*180/pi,'%.1f') '^{\circ}, A1 = ' num2str(a_1*180/pi,'%.1f') '^{\circ}'],'FontWeight','normal','FontSize',axes_size*0.6);
elseif strcmp(authors{1},'GU') && length(authors) == 1
    title(['GUD ' num2str(GUD) ': M = ' num2str(M,'%.3f') ', k = ' num2str(k,'%.2f') ', A0 = ' num2str(a_0*180/pi,'%.1f') '^{\circ}, A1 = ' num2str(a_1*180/pi,'%.1f') '^{\circ}'],'FontWeight','normal','FontSize',axes_size*0.6);
end
grid on

% c_c
figure7 = figure('InvertHardcopy','off','Color',[1 1 1]); figure7.PaperSize = [6 4.5];
axes7 = axes('Parent',figure7,'FontSize',axes_size,'FontName','times new roman');
axes7.XLim = xlim_vec;
hold(axes7,'on');
if length(authors) == 1 
    if ~isnan(alpha_interp(1)) && strcmp(authors{1},'McAlister et al. (1982)')
        plot7 = plot(alpha(iti:end)*180/pi,c_c(iti:end),'k-',alpha_BLS(iti_BLS:end)*180/pi,cc_BLS(iti_BLS:end),'k--',alpha_interp*180/pi,cc_interp,'b-','LineWidth',lw,'Parent',axes7);
        legend('Present model','Sheng et al. model',[authors{1}]);
    elseif isnan(data{15}(1)) || isnan(authors{1}(1))
        plot7 = plot(alpha(iti:end)*180/pi,c_c(iti:end),'k-',alpha_BLS(iti_BLS:end)*180/pi,cc_BLS(iti_BLS:end),'k--','LineWidth',lw,'Parent',axes7);
        legend('Present model','Sheng et al. model');
    else
        plot7 = plot(alpha(iti:end)*180/pi,c_c(iti:end),'k-',alpha_BLS(iti_BLS:end)*180/pi,cc_BLS(iti_BLS:end),'k--',data{15,1},data{16,1},'ko','LineWidth',lw,'MarkerSize',ms,'Parent',axes7);
        legend('Present model','Sheng et al. model',authors{1});
    end
else
    plot7 = plot(alpha(iti:end)*180/pi,c_c(iti:end),'k-',data{15,1},data{16,1},'ko',data{25,1},data{26,1},'k--','LineWidth',lw,'MarkerSize',ms,'Parent',axes7);
    if isnan(data{15}(1))
        legend('Present model');
    else
        legend('Present model',authors{1},authors{2});
    end
end
xlabel('\alpha [deg]','FontWeight','normal','FontSize',axes_size);
ylabel('Chordwise coefficient','FontWeight','normal','FontSize',axes_size);
if strcmp(authors{1},'McAlister et al. (1982)') && length(authors) == 1
    title(['Frame ' num2str(frame) ': M = ' num2str(M) ', k = ' num2str(k) ', A0 = ' num2str(a_0*180/pi,'%.1f') '^{\circ}, A1 = ' num2str(a_1*180/pi,'%.1f') '^{\circ}'],'FontWeight','normal','FontSize',axes_size*0.6);
elseif strcmp(authors{1},'GU') && length(authors) == 1
    title(['GUD ' num2str(GUD) ': M = ' num2str(M,'%.3f') ', k = ' num2str(k,'%.2f') ', A0 = ' num2str(a_0*180/pi,'%.1f') '^{\circ}, A1 = ' num2str(a_1*180/pi,'%.1f') '^{\circ}'],'FontWeight','normal','FontSize',axes_size*0.6);
end
grid on
set(legend,'Location','best','FontSize',axes_size*0.6);

% c_l x time
figure4 = figure('InvertHardcopy','off','Color',[1 1 1]); figure4.PaperSize = [6 4.5];
axes4 = axes('Parent',figure4,'FontSize',axes_size,'FontName','times new roman','XTick',[-90 0 90 180 270]);
axes4.XLim = [-90 270];
hold(axes4,'on'); 
if isnan(data{7}(1)) || isnan(authors{1}(1))
    plot4 = plot((t(itt:ittf)-t(itt))/t_cycle*360-90,c_l(itt:ittf),'k-',(t_BLS(itt_BLS:ittf_BLS)-t_BLS(itt_BLS))/t_cycle*360-90,cl_BLS(itt_BLS:ittf_BLS),'k--','LineWidth',lw,'Parent',axes4);
    legend('Present model','Sheng et al. model');
else
    plot4 = plot((t(itt:ittf)-t(itt))/t_cycle*360-90,c_l(itt:ittf),'k-',(t_BLS(itt_BLS:ittf_BLS)-t_BLS(itt_BLS))/t_cycle*360-90,cl_BLS(itt_BLS:ittf_BLS),'k--',data{7,1},data{8,1},'ko','LineWidth',lw,'MarkerSize',ms,'Parent',axes4);
    legend('Present model','Sheng et al. model',authors{1});
end
xlabel('\omega t [deg]','FontWeight','normal','FontSize',axes_size);
ylabel('Lift coefficient','FontWeight','normal','FontSize',axes_size);
if strcmp(authors{1},'McAlister et al. (1982)') && length(authors) == 1
    title(['Frame ' num2str(frame) ': M = ' num2str(M) ', k = ' num2str(k) ', A0 = ' num2str(a_0*180/pi,'%.1f') '^{\circ}, A1 = ' num2str(a_1*180/pi,'%.1f') '^{\circ}'],'FontWeight','normal','FontSize',axes_size*0.6);
elseif strcmp(authors{1},'GU') && length(authors) == 1
    title(['GUD ' num2str(GUD) ': M = ' num2str(M,'%.3f') ', k = ' num2str(k,'%.2f') ', A0 = ' num2str(a_0*180/pi,'%.1f') '^{\circ}, A1 = ' num2str(a_1*180/pi,'%.1f') '^{\circ}'],'FontWeight','normal','FontSize',axes_size*0.6);
end
set(legend,'Location','best','FontSize',axes_size*0.6);
grid on

% c_m x time
figure5 = figure('InvertHardcopy','off','Color',[1 1 1]); figure5.PaperSize = [6 4.5];
axes5 = axes('Parent',figure5,'FontSize',axes_size,'FontName','times new roman','XTick',[-90 0 90 180 270]);
axes5.XLim = [-90 270];
hold(axes5,'on'); 
if isnan(data{9}(1)) || isnan(authors{1}(1))
    plot5 = plot((t(itt:ittf)-t(itt))/t_cycle*360-90,c_m(itt:ittf),'k-',(t_BLS(itt_BLS:ittf_BLS)-t_BLS(itt_BLS))/t_cycle*360-90,cm_BLS(itt_BLS:ittf_BLS),'k--','LineWidth',lw,'Parent',axes5);
    legend('Present model','Sheng et al. model');
else
    plot5 = plot((t(itt:ittf)-t(itt))/t_cycle*360-90,c_m(itt:ittf),'k-',(t_BLS(itt_BLS:ittf_BLS)-t_BLS(itt_BLS))/t_cycle*360-90,cm_BLS(itt_BLS:ittf_BLS),'k--',data{9,1},data{10,1},'ko','LineWidth',lw,'MarkerSize',ms,'Parent',axes5);
    legend('Present model','Sheng et al. model',authors{1});
end
xlabel('\omega t [deg]','FontWeight','normal','FontSize',axes_size);
ylabel('Moment coefficient','FontWeight','normal','FontSize',axes_size);
if strcmp(authors{1},'McAlister et al. (1982)') && length(authors) == 1
    title(['Frame ' num2str(frame) ': M = ' num2str(M) ', k = ' num2str(k) ', A0 = ' num2str(a_0*180/pi,'%.1f') '^{\circ}, A1 = ' num2str(a_1*180/pi,'%.1f') '^{\circ}'],'FontWeight','normal','FontSize',axes_size*0.6);
elseif strcmp(authors{1},'GU') && length(authors) == 1
    title(['GUD ' num2str(GUD) ': M = ' num2str(M,'%.3f') ', k = ' num2str(k,'%.2f') ', A0 = ' num2str(a_0*180/pi,'%.1f') '^{\circ}, A1 = ' num2str(a_1*180/pi,'%.1f') '^{\circ}'],'FontWeight','normal','FontSize',axes_size*0.6);
end
set(legend,'Location','best','FontSize',axes_size*0.6);
grid on

% c_d x time
figure6 = figure('InvertHardcopy','off','Color',[1 1 1]); figure6.PaperSize = [6 4.5];
axes6 = axes('Parent',figure6,'FontSize',axes_size,'FontName','times new roman','XTick',[-90 0 90 180 270]);
axes6.XLim = [-90 270];
hold(axes6,'on'); 
if isnan(data{11}(1)) || isnan(authors{1}(1))
    plot6 = plot((t(itt:ittf)-t(itt))/t_cycle*360-90,c_d(itt:ittf),'k-',(t_BLS(itt_BLS:ittf_BLS)-t_BLS(itt_BLS))/t_cycle*360-90,cd_BLS(itt_BLS:ittf_BLS),'k--','LineWidth',lw,'Parent',axes6);
    legend('Present model','Sheng et al. model');
else
    plot6 = plot((t(itt:ittf)-t(itt))/t_cycle*360-90,c_d(itt:ittf),'k-',(t_BLS(itt_BLS:ittf_BLS)-t_BLS(itt_BLS))/t_cycle*360-90,cd_BLS(itt_BLS:ittf_BLS),'k--',data{11,1},data{12,1},'ko','LineWidth',lw,'MarkerSize',ms,'Parent',axes6);
    legend('Present model','Sheng et al. model',authors{1});
end
xlabel('\omega t [deg]','FontWeight','normal','FontSize',axes_size);
ylabel('Drag coefficient','FontWeight','normal','FontSize',axes_size);
if strcmp(authors{1},'McAlister et al. (1982)') && length(authors) == 1
    title(['Frame ' num2str(frame) ': M = ' num2str(M) ', k = ' num2str(k) ', A0 = ' num2str(a_0*180/pi,'%.1f') '^{\circ}, A1 = ' num2str(a_1*180/pi,'%.1f') '^{\circ}'],'FontWeight','normal','FontSize',axes_size*0.6);
elseif strcmp(authors{1},'GU') && length(authors) == 1
    title(['GUD ' num2str(GUD) ': M = ' num2str(M,'%.3f') ', k = ' num2str(k,'%.2f') ', A0 = ' num2str(a_0*180/pi,'%.1f') '^{\circ}, A1 = ' num2str(a_1*180/pi,'%.1f') '^{\circ}'],'FontWeight','normal','FontSize',axes_size*0.6);
end
set(legend,'Location','best','FontSize',axes_size*0.6);
grid on

% c_n x time
figure10 = figure('InvertHardcopy','off','Color',[1 1 1]); figure10.PaperSize = [6 4.5];
axes10 = axes('Parent',figure10,'FontSize',axes_size,'FontName','times new roman','XTick',[-90 0 90 180 270]);
axes10.XLim = [-90 270];
hold(axes10,'on');
if ~isnan(cn_interp(1)) && strcmp(authors{1},'McAlister et al. (1982)')
    plot10 = plot((t(itt:ittf)-t(itt))/t_cycle*360-90,c_n(itt:ittf),'k-',(t_BLS(itt_BLS:ittf_BLS)-t_BLS(itt_BLS))/t_cycle*360-90,cn_BLS(itt_BLS:ittf_BLS),'k--',time_interp,cn_interp,'b-','LineWidth',lw,'Parent',axes10);
    legend('Present model','Sheng et al. model',[authors{1}]);
elseif ~isnan(cn_interp(1)) && strcmp(authors{1},'GU')
    plot10 = plot((t(itt:ittf)-t(itt))/t_cycle*360-90,c_n(itt:ittf),'k-',(t_BLS(itt_BLS:ittf_BLS)-t_BLS(itt_BLS))/t_cycle*360-90,cn_BLS(itt_BLS:ittf_BLS),'k--',data{27,1},data{28,1},'ko','LineWidth',lw,'Parent',axes10);
    legend('Present model','Sheng et al. model',[authors{1}]);
else
    plot10 = plot((t(itt:ittf)-t(itt))/t_cycle*360-90,c_n(itt:ittf),'k-',(t_BLS(itt_BLS:ittf_BLS)-t_BLS(itt_BLS))/t_cycle*360-90,cn_BLS(itt_BLS:ittf_BLS),'k--','LineWidth',lw,'Parent',axes10);
    legend('Present model','Sheng et al. model');
end
xlabel('\omega t [deg]','FontWeight','normal','FontSize',axes_size);
ylabel('Normal coefficient','FontWeight','normal','FontSize',axes_size);
if strcmp(authors{1},'McAlister et al. (1982)') && length(authors) == 1
    title(['Frame ' num2str(frame) ': M = ' num2str(M) ', k = ' num2str(k) ', A0 = ' num2str(a_0*180/pi,'%.1f') '^{\circ}, A1 = ' num2str(a_1*180/pi,'%.1f') '^{\circ}'],'FontWeight','normal','FontSize',axes_size*0.6);
elseif strcmp(authors{1},'GU') && length(authors) == 1
    title(['GUD ' num2str(GUD) ': M = ' num2str(M,'%.3f') ', k = ' num2str(k,'%.2f') ', A0 = ' num2str(a_0*180/pi,'%.1f') '^{\circ}, A1 = ' num2str(a_1*180/pi,'%.1f') '^{\circ}'],'FontWeight','normal','FontSize',axes_size*0.6);
end
set(legend,'Location','best','FontSize',axes_size*0.6);
grid on

% c_c x time
figure9 = figure('InvertHardcopy','off','Color',[1 1 1]); figure9.PaperSize = [6 4.5];
axes9 = axes('Parent',figure9,'FontSize',axes_size,'FontName','times new roman','XTick',[-90 0 90 180 270]);
axes9.XLim = [-90 270];
hold(axes9,'on'); 
if ~isnan(cc_interp(1)) && strcmp(authors{1},'McAlister et al. (1982)')
    plot9 = plot((t(itt:ittf)-t(itt))/t_cycle*360-90,c_c(itt:ittf),'k-',(t_BLS(itt_BLS:ittf_BLS)-t_BLS(itt_BLS))/t_cycle*360-90,cc_BLS(itt_BLS:ittf_BLS),'k--',time_interp,cc_interp,'b-','LineWidth',lw,'Parent',axes9);
    legend('Present model','Sheng et al. model',[authors{1}]);
elseif ~isnan(cn_interp(1)) && strcmp(authors{1},'GU')
    plot9 = plot((t(itt:ittf)-t(itt))/t_cycle*360-90,c_c(itt:ittf),'k-',(t_BLS(itt_BLS:ittf_BLS)-t_BLS(itt_BLS))/t_cycle*360-90,cc_BLS(itt_BLS:ittf_BLS),'k--',data{29,1},data{30,1},'ko','LineWidth',lw,'Parent',axes9);
    legend('Present model','Sheng et al. model',[authors{1}]);
else
    plot9 = plot((t(itt:ittf)-t(itt))/t_cycle*360-90,c_c(itt:ittf),'k-',(t_BLS(itt_BLS:ittf_BLS)-t_BLS(itt_BLS))/t_cycle*360-90,cc_BLS(itt_BLS:ittf_BLS),'k--','LineWidth',lw,'Parent',axes9);
    legend('Present model','Sheng et al. model');
end
xlabel('\omega t [deg]','FontWeight','normal','FontSize',axes_size);
ylabel('Chordwise coefficient','FontWeight','normal','FontSize',axes_size);
if strcmp(authors{1},'McAlister et al. (1982)') && length(authors) == 1
    title(['Frame ' num2str(frame) ': M = ' num2str(M) ', k = ' num2str(k) ', A0 = ' num2str(a_0*180/pi,'%.1f') '^{\circ}, A1 = ' num2str(a_1*180/pi,'%.1f') '^{\circ}'],'FontWeight','normal','FontSize',axes_size*0.6);
elseif strcmp(authors{1},'GU') && length(authors) == 1
    title(['GUD ' num2str(GUD) ': M = ' num2str(M,'%.3f') ', k = ' num2str(k,'%.2f') ', A0 = ' num2str(a_0*180/pi,'%.1f') '^{\circ}, A1 = ' num2str(a_1*180/pi,'%.1f') '^{\circ}'],'FontWeight','normal','FontSize',axes_size*0.6);
end
set(legend,'Location','best','FontSize',axes_size*0.6);
grid on

% resave_figures;