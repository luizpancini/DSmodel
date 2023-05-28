clc
clear 
% close all
addpath('../functions')

M = 0.3;
% Only stalling frames
if M == 0.3
    frame_list = [7104; 7019; 7021; 7101; 7023; 7112; 7222; 7302; 7113; 7300; 7114; 10202; 10203; 10204; 10207; 10208; 10211; 10212; 10303; 7117; 7118; 7119; 7120; 7121; 7200; 7202; 7205; 7207; 7305; 9221; 9223; 9302; 9307; 10104; 10105; 10108; 9210; 9213; 9214; 10022; 10113; 10114; 10117; 10118; 10120; 10123; 14218; 14219];
    t_min_cn =   [90;   90;   150;  90;   90;   90;   150;  140;  150;  150;  140;  120;   140;   150;   150;   150;   150;   150;   150;   150;  150;  150;  150;  150;  150;  150;  150;  150;  150;  140;  150;  150;  150;  150;   150;   140;   150;  150;  130;  150;   150;   150;   150;   150;   150;   150;   150;   150];
elseif M == 0.185
    frame_list = [8220; 8222; 8306; 9022; 9110; 9112];
end

j = 0;
N = length(frame_list);
for i=1:N
    frame = frame_list(i);
    [b,a_inf,airfoil,M,alpha_0,delta_alpha,k,alpha_exp_cl,cl_exp,alpha_exp_cm,cm_exp,alpha_exp_cd,cd_exp,alpha_exp_cn,cn_exp,alpha_exp_cc,cc_exp,time_exp_cl,clt_exp,time_exp_cm,cmt_exp,time_exp_cd,cdt_exp,authors] = load_frame(frame);
    data = {alpha_exp_cl; cl_exp; alpha_exp_cm; cm_exp; alpha_exp_cd; cd_exp; time_exp_cl; clt_exp; time_exp_cm; cmt_exp; time_exp_cd; cdt_exp; alpha_exp_cn; cn_exp; alpha_exp_cc; cc_exp};
    dt = 2e-4;                    % Maximum time step
    t_cycle = 2*pi*b/(M*a_inf*k); % Time of one cycle
    [time_interp, alpha_interp, cn_interp, cm_interp, cc_interp] = interp_coefs(data,alpha_0,delta_alpha,round(t_cycle/dt));
    [alpha_re_maxcm_exp,r_re_maxcm_exp,alpha_re_mincn_exp,r_re_mincn_exp,alpha_re_maxcm,r_re_maxcm,alpha_re_mincn,r_re_mincn,alpha_re_a1m,r_re_a1m] = find_re_data(time_interp,cn_interp,cm_interp,alpha_0,delta_alpha,k,M*a_inf,b,M,frame,t_min_cn,i);
    alpha_re_maxcm_exp_vec(i) = alpha_re_maxcm_exp;
    r_re_exp_vec(i) = r_re_maxcm_exp;
    alpha_re_mincn_exp_vec(i) = alpha_re_mincn_exp;
    r_re_mincn_exp_vec(i) = r_re_mincn_exp;
    alpha_re_maxcm_vec(i) = alpha_re_maxcm;
    r_re_maxcm_vec(i) = r_re_maxcm;
    alpha_re_mincn_vec(i) = alpha_re_mincn;
    r_re_mincn_vec(i) = r_re_mincn;
    alpha_re_a1m_vec(i) = alpha_re_a1m;
    r_re_a1m_vec(i) = r_re_a1m;
    k_vec(i) = k;
end

% Sheng et al.'s criterion
% For M = 0.3
alpha_min0 = 18.5*pi/180;
T_r = 7.0;
r_f = 0.025;
x_dashed = [0 r_f];
y_dashed = [alpha_min0 alpha_min0-T_r*r_f];

%% Plots
set(0,'DefaultTextInterpreter','tex')
set(0,'DefaultLegendInterpreter','tex')
axes_size = 20;
lw = 0.75;
ms = 6;

figure1 = figure('InvertHardcopy','off','Color',[1 1 1]);
axes1 = axes('Parent',figure1,'FontSize',axes_size,'FontName','times new roman');
hold(axes1,'on'); 
p1 = plot(-r_re_mincn_exp_vec,alpha_re_mincn_exp_vec*180/pi,'ks','MarkerFaceColor','k','MarkerSize',ms,'Parent',axes1);
p2 = plot(-r_re_mincn_vec,alpha_re_mincn_vec*180/pi,'bo','MarkerSize',ms,'Parent',axes1);
% p3 = plot(-x_dashed,y_dashed*180/pi,'b--','LineWidth',lw,'Parent',axes1);
xlabel('r','FontWeight','normal','FontSize',axes_size);
ylabel('\alpha_{re} [deg]','FontWeight','normal','FontSize',axes_size);
ylim([5 20]);
grid on
lgd1 = legend([p1 p2],'Experimental','Model');
% lgd1 = legend([p1 p2 p3],'Experimental','Model','Sheng et al. fit');
set(lgd1,'Location','northwest','FontSize',14,'Box','on');
figure1.PaperSize = [6 4.5];

figure2 = figure('InvertHardcopy','off','Color',[1 1 1]); figure2.PaperSize = [6 4.5];
axes2 = axes('Parent',figure2,'FontSize',axes_size,'FontName','times new roman'); axes2.XDir = 'reverse';
hold(axes2,'on'); 
c = jet(length(frame_list));
p1 = plot(nan,nan,'kx','MarkerFaceColor','none','Parent',axes2);
p2 = plot(nan,nan,'ko','MarkerFaceColor','none','Parent',axes2);
for i=1:length(frame_list)
    plot(-r_re_mincn_vec(i),alpha_re_mincn_vec(i)*180/pi,'o','MarkerFaceColor',c(i,:),'MarkerEdgeColor',c(i,:),'MarkerSize',ms,'Parent',axes2);
    plot(-r_re_mincn_exp_vec(i),alpha_re_mincn_exp_vec(i)*180/pi,'x','MarkerFaceColor',c(i,:),'MarkerEdgeColor',c(i,:),'MarkerSize',1.5*ms,'Parent',axes2);
end
xlabel('r','FontWeight','normal','FontSize',axes_size);
ylabel('\alpha_{re} [deg]','FontWeight','normal','FontSize',axes_size);
ylim([5 20]);
grid on
lgd1 = legend([p1 p2],'Experimental','Model');
set(lgd1,'Location','northeast','FontSize',14,'Box','on','AutoUpdate','off');
scatter(nan,nan,nan,nan)
colormap(jet(8));
cb = colorbar(); cb.Limits = [0 1]; cb.Ticks = [0,1]; cb.TickLabels = {'Light stall','Deep stall'}; cb.Location = 'northoutside';

figure3 = figure('InvertHardcopy','off','Color',[1 1 1]); figure3.PaperSize = [6 4.5];
axes3 = axes('Parent',figure3,'FontSize',axes_size,'FontName','times new roman'); 
hold(axes3,'on'); 
c = jet(length(frame_list));
p1 = plot(nan,nan,'kx','MarkerFaceColor','none','Parent',axes3);
p2 = plot(nan,nan,'ko','MarkerFaceColor','none','Parent',axes3);
for i=1:length(frame_list)
    plot(k_vec(i),alpha_re_mincn_vec(i)*180/pi,'o','MarkerFaceColor',c(i,:),'MarkerEdgeColor',c(i,:),'MarkerSize',ms,'Parent',axes3);
    plot(k_vec(i),alpha_re_mincn_exp_vec(i)*180/pi,'x','MarkerFaceColor',c(i,:),'MarkerEdgeColor',c(i,:),'MarkerSize',1.5*ms,'Parent',axes3);
end
xlabel('k','FontWeight','normal','FontSize',axes_size);
ylabel('\alpha_{re} [deg]','FontWeight','normal','FontSize',axes_size);
xlim([0 inf]); ylim([5 20]);
grid on
lgd1 = legend([p1 p2],'Experimental','Model');
set(lgd1,'Location','northeast','FontSize',14,'Box','on','AutoUpdate','off');
scatter(nan,nan,nan,nan)
colormap(jet(8));
cb = colorbar(); cb.Limits = [0 1]; cb.Ticks = [0,1]; cb.TickLabels = {'Light stall','Deep stall'}; cb.Location = 'northoutside';

% figure3 = figure('InvertHardcopy','off','Color',[1 1 1]);
% axes3 = axes('Parent',figure3,'FontSize',axes_size,'FontName','times new roman');
% hold(axes3,'on'); 
% p1 = plot(alpha_re_a1m_vec*180/pi,alpha_re_mincn_exp_vec*180/pi,'ko',[6 16],[6 16],'b--',[8 14],[6.5 12.5],'k--','MarkerSize',ms,'Parent',axes3);
% xlabel('\alpha_{re} [deg] - \alpha^{\prime} = \alpha_{1_m} criterion','FontWeight','normal','FontSize',axes_size);
% ylabel('\alpha_{re} [deg] - min c_n criterion','FontWeight','normal','FontSize',axes_size);
% grid on
% figure3.PaperSize = [6 4.5];
