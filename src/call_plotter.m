function h = call_plotter(INPUTS,case_now,data,params,outputs,interp,varargin)

% Handle inputs
if isempty(varargin)
    plot_mod_vars = true;
else
    plot_mod_vars = varargin{1};
end

% Initialize plot handles structure
h = struct;

% Unpack 
authors = data.authors;
model = INPUTS.model;
source = INPUTS.source; if isfield(data,'source'),source = data.source; end
base_name = INPUTS.base_name;
t = outputs.t;
alpha = outputs.alpha;
c_n = outputs.c_n;
c_m = outputs.c_m;
c_c = outputs.c_c;
c_l = outputs.c_l;
c_d = outputs.c_d;
iti = interp.iti;
time_ref_interp = interp.time_ref;
alpha_ref_interp = interp.alpha_ref;
cn_ref_interp = interp.cn_ref;
cm_ref_interp = interp.cm_ref;
cc_ref_interp = interp.cc_ref;
cl_ref_interp = interp.cl_ref;
cd_ref_interp = interp.cd_ref;

%% Set legend for current model
switch model
    case {"BL","BLT","BLgust","BLhargen","BLTgust","BLTflap","BLtvfs","BLTtvfs","BLThargen"}
        model_name = "Present model";
    case {"BLS","BLSLR"}
        model_name = "Sheng et al. model";
    case {"BLO","BLOgust"}
        model_name = "Original model";   
    case "BLG"
        model_name = "Gupta & Leishman model";        
end

%% Create figure and initialize plot options
[h,plot_opt] = init_plotter(h,base_name,case_now,source,params,1);

%% Coefficients x pitch angle
if ~contains(model,["gust","flap","tvfs"])
    h_cla = alpha_plotter(1,h.fig,'cl','$c_l$',alpha,c_l,alpha_ref_interp,cl_ref_interp,data.alpha_exp_cl,data.cl_exp,data.alpha_mod_cl,data.cl_mod,data.alpha_cycles,data.cl_cycles,authors,iti,source,h.tabgp,model_name,plot_opt);
    h_cma = alpha_plotter(1,h.fig,'cm','$c_m$',alpha,c_m,alpha_ref_interp,cm_ref_interp,data.alpha_exp_cm,data.cm_exp,data.alpha_mod_cm,data.cm_mod,data.alpha_cycles,data.cm_cycles,authors,iti,source,h.tabgp,model_name,plot_opt);
    h_cda = alpha_plotter(1,h.fig,'cd','$c_d$',alpha,c_d,alpha_ref_interp,cd_ref_interp,data.alpha_exp_cd,data.cd_exp,data.alpha_mod_cd,data.cd_mod,data.alpha_cycles,data.cd_cycles,authors,iti,source,h.tabgp,model_name,plot_opt);
    h_cna = alpha_plotter(1,h.fig,'cn','$c_n$',alpha,c_n,alpha_ref_interp,cn_ref_interp,data.alpha_exp_cn,data.cn_exp,data.alpha_mod_cn,data.cn_mod,data.alpha_cycles,data.cn_cycles,authors,iti,source,h.tabgp,model_name,plot_opt);
    h_cca = alpha_plotter(1,h.fig,'cc','$c_c$',alpha,c_c,alpha_ref_interp,cc_ref_interp,data.alpha_exp_cc,data.cc_exp,data.alpha_mod_cc,data.cc_mod,data.alpha_cycles,data.cc_cycles,authors,iti,source,h.tabgp,model_name,plot_opt);
    % Add to plots handle structure
    h = variables2struct(h,h_cla,h_cma,h_cda,h_cna,h_cca);
end

%% Coefficients x plunge-induced AoA
if isfield(outputs,'alpha_plunge')
    alpha_plunge = outputs.alpha_plunge;
    h_clap = alpha_plunge_plotter(1,h.fig,'cl','$c_l$',alpha_plunge,c_l,data.alpha_exp_cl,data.cl_exp,data.alpha_mod_cl,data.cl_mod,authors,iti,h.tabgp,model_name,plot_opt);
    h_cmap = alpha_plunge_plotter(1,h.fig,'cm','$c_m$',alpha_plunge,c_m,data.alpha_exp_cm,data.cm_exp,data.alpha_mod_cm,data.cm_mod,authors,iti,h.tabgp,model_name,plot_opt);
    h_cdap = alpha_plunge_plotter(1,h.fig,'cd','$c_d$',alpha_plunge,c_d,data.alpha_exp_cd,data.cd_exp,data.alpha_mod_cd,data.cd_mod,authors,iti,h.tabgp,model_name,plot_opt);
    h_cnap = alpha_plunge_plotter(1,h.fig,'cn','$c_n$',alpha_plunge,c_n,data.alpha_exp_cn,data.cn_exp,data.alpha_mod_cn,data.cn_mod,authors,iti,h.tabgp,model_name,plot_opt);
    h_ccap = alpha_plunge_plotter(1,h.fig,'cc','$c_c$',alpha_plunge,c_c,data.alpha_exp_cc,data.cc_exp,data.alpha_mod_cc,data.cc_mod,authors,iti,h.tabgp,model_name,plot_opt);
   % Add to plots handle structure
    h = variables2struct(h,h_clap,h_cmap,h_cdap,h_cnap,h_ccap);
end

%% Coefficients x flap angle
if contains(model,["flap","gen"])
    delta = outputs.delta;
    c_h = outputs.c_h;
    h_cld = delta_plotter('Lift coefficient','cl',delta,c_l,data.delta_exp_cl,data.cl_exp,data.delta_mod_cl,data.cl_mod,authors,iti,h.tabgp,model_name,plot_opt);
    h_cmd = delta_plotter('Moment coefficient','cm',delta,c_m,data.delta_exp_cm,data.cm_exp,data.delta_mod_cm,data.cm_mod,authors,iti,h.tabgp,model_name,plot_opt);
    h_cdd = delta_plotter('Drag coefficient','cd',delta,c_d,data.delta_exp_cd,data.cd_exp,data.delta_mod_cd,data.cd_mod,authors,iti,h.tabgp,model_name,plot_opt);
    h_cnd = delta_plotter('Normal coefficient','cn',delta,c_n,data.delta_exp_cn,data.cn_exp,data.delta_mod_cn,data.cn_mod,authors,iti,h.tabgp,model_name,plot_opt);
    h_ccd = delta_plotter('Chordwise coefficient','cc',delta,c_c,data.delta_exp_cc,data.cc_exp,data.delta_mod_cc,data.cc_mod,authors,iti,h.tabgp,model_name,plot_opt);
    h_chd = delta_plotter('Hinge coefficient','ch',delta,c_h,data.delta_exp_ch,data.ch_exp,data.delta_mod_ch,data.ch_mod,authors,iti,h.tabgp,model_name,plot_opt);
    % Add to plots handle structure
    h = variables2struct(h,h_cld,h_cmd,h_cdd,h_cnd,h_ccd,h_chd);
end

%% Coefficients x time
h_clt = time_plotter(1,h.fig,'cl x t','$c_l$',t,c_l,time_ref_interp,cl_ref_interp,data.time_exp_cl,data.clt_exp,data.time_mod_cl,data.clt_mod,data.time_cycles,data.cl_cycles,model,outputs,params,interp,authors,source,h.tabgp,model_name,plot_opt);
h_cmt = time_plotter(1,h.fig,'cm x t','$c_m$',t,c_m,time_ref_interp,cm_ref_interp,data.time_exp_cm,data.cmt_exp,data.time_mod_cm,data.cmt_mod,data.time_cycles,data.cm_cycles,model,outputs,params,interp,authors,source,h.tabgp,model_name,plot_opt);
h_cdt = time_plotter(1,h.fig,'cd x t','$c_d$',t,c_d,time_ref_interp,cd_ref_interp,data.time_exp_cd,data.cdt_exp,data.time_mod_cd,data.cdt_mod,data.time_cycles,data.cd_cycles,model,outputs,params,interp,authors,source,h.tabgp,model_name,plot_opt);
h_cnt = time_plotter(1,h.fig,'cn x t','$c_n$',t,c_n,time_ref_interp,cn_ref_interp,data.time_exp_cn,data.cnt_exp,data.time_mod_cn,data.cnt_mod,data.time_cycles,data.cn_cycles,model,outputs,params,interp,authors,source,h.tabgp,model_name,plot_opt);
h_cct = time_plotter(1,h.fig,'cc x t','$c_c$',t,c_c,time_ref_interp,cc_ref_interp,data.time_exp_cc,data.cct_exp,data.time_mod_cc,data.cct_mod,data.time_cycles,data.cc_cycles,model,outputs,params,interp,authors,source,h.tabgp,model_name,plot_opt);
% Add to plots handle structure
h = variables2struct(h,h_clt,h_cmt,h_cdt,h_cnt,h_cct);
    
%% Separation points, angle offsets and time delay constants
if ismember(model,{'BL','BLT'}) && plot_mod_vars
    if ~isnan(alpha_ref_interp)
        f_da_Tf_plotter(interp,outputs,source,model,h.tabgp,plot_opt,params);
    end
end

end