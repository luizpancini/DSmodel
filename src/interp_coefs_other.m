function [time,alpha,cn_interp,cm_interp,cc_interp,cl_interp,cd_interp,interp] = interp_coefs_other(case_now,data,a_0,a_1,t,t_cycle,interp,N_interp,interp_method)

if ismember(case_now,11:31)
    % Unpack
    alpha_exp_cn = data.alpha_exp_cn;
    alpha_exp_cm = data.alpha_exp_cm;
    cn_exp = data.cn_exp;
    cm_exp = data.cm_exp;
    
    % Find time shift and new indices
    alpha_init_n = alpha_exp_cn(1)*pi/180;                  % Pitch angle at which the experiment begins
    t_init_n = asind((alpha_init_n-a_0)/a_1);               % Cycle angle [deg] at which the data starts
    delta_t_n = -t_init_n/360;                              % According time shift (fraction of t_cycle)
    itt_n = find(t>t(end)-(1+delta_t_n)*t_cycle,1,'first');
    ittf_n = find(t>t(end)-(delta_t_n)*t_cycle,1,'first');
    alpha_init_m = alpha_exp_cm(1)*pi/180;                  % Pitch angle at which the experiment begins
    t_init_m = asind((alpha_init_m-a_0)/a_1);               % Cycle angle [deg] at which the data starts
    delta_t_m = -t_init_m/360;                              % According time shift (fraction of t_cycle)
    itt_m = find(t>t(end)-(1+delta_t_m)*t_cycle,1,'first');
    ittf_m = find(t>t(end)-(delta_t_m)*t_cycle,1,'first');
    
    % Update interpolation data structure
    interp.itt = min([itt_n,itt_m]);
    interp.ittf = max([ittf_n,ittf_m]);
    interp.range = interp.itt:interp.ittf;
    
    % Cycle angle (time) and AoA vectors
    time = linspace(-90,270,N_interp);
    alpha = 180/pi*(a_0+a_1*sind(time));
    [~,ind_alpha_exp_cn_max] = max(alpha_exp_cn);
    [~,ind_alpha_exp_cm_max] = max(alpha_exp_cm);
    a_0_n = (max(alpha_exp_cn)+min(alpha_exp_cn))/2;
    a_1_n = (max(alpha_exp_cn)-min(alpha_exp_cn))/2;
    a_0_m = (max(alpha_exp_cm)+min(alpha_exp_cm))/2;
    a_1_m = (max(alpha_exp_cm)-min(alpha_exp_cm))/2;
    time_exp_cn = real(asind((alpha_exp_cn-a_0_n)/a_1_n));
    time_exp_cn(ind_alpha_exp_cn_max+1:end) = 180-time_exp_cn(ind_alpha_exp_cn_max+1:end);
    time_exp_cm = real(asind((alpha_exp_cm-a_0_m)/a_1_m));
    time_exp_cm(ind_alpha_exp_cm_max+1:end) = 180-time_exp_cm(ind_alpha_exp_cm_max+1:end);
    
    % Interpolate coefficients
    cn_interp = interp1(time_exp_cn,cn_exp,time,interp_method);
    cm_interp = interp1(time_exp_cm,cm_exp,time,interp_method);
    cc_interp = nan*cn_interp; cl_interp = nan*cn_interp; cd_interp = nan*cn_interp;
        
else
    time = nan;
    alpha = nan;
    cn_interp = nan;
    cm_interp = nan;
    cc_interp = nan;
    cl_interp = nan;
    cd_interp = nan;
end

end