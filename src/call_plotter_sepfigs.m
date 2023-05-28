function handles = call_plotter_sepfigs(INPUTS,case_now,data,params,outputs,interp,save_folder,figures_extension)

%% Handle inputs
% Set figures extension
figures_extension = string(figures_extension);

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

% Set filepaths
baseFileName_cn_a = base_name + num2str(case_now) + "_cn_a" + figures_extension;
baseFileName_cm_a = base_name + num2str(case_now) + "_cm_a" + figures_extension;
baseFileName_cc_a = base_name + num2str(case_now) + "_cc_a" + figures_extension;
baseFileName_cl_a = base_name + num2str(case_now) + "_cl_a" + figures_extension;
baseFileName_cd_a = base_name + num2str(case_now) + "_cd_a" + figures_extension;
baseFileName_cn_t = base_name + num2str(case_now) + "_cn_t" + figures_extension;
baseFileName_cm_t = base_name + num2str(case_now) + "_cm_t" + figures_extension;
baseFileName_cc_t = base_name + num2str(case_now) + "_cc_t" + figures_extension;
baseFileName_cl_t = base_name + num2str(case_now) + "_cl_t" + figures_extension;
baseFileName_cd_t = base_name + num2str(case_now) + "_cd_t" + figures_extension;
filepath_cn_a = fullfile(save_folder,baseFileName_cn_a);
filepath_cm_a = fullfile(save_folder,baseFileName_cm_a);
filepath_cc_a = fullfile(save_folder,baseFileName_cc_a);
filepath_cl_a = fullfile(save_folder,baseFileName_cl_a);
filepath_cd_a = fullfile(save_folder,baseFileName_cd_a);
filepath_cn_t = fullfile(save_folder,baseFileName_cn_t);
filepath_cm_t = fullfile(save_folder,baseFileName_cm_t);
filepath_cc_t = fullfile(save_folder,baseFileName_cc_t);
filepath_cl_t = fullfile(save_folder,baseFileName_cl_t);
filepath_cd_t = fullfile(save_folder,baseFileName_cd_t);

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
% Initialize plot handles structure array
for j=1:16
    [handles(j),plot_opt] = init_plotter(struct,base_name,case_now,source,params,0);
end

%% Coefficients x pitch angle
if ~contains(model,["gust","flap","tvfs"])
    [~,handles(1).fig] = alpha_plotter(0,handles(1).fig,'cl','$c_l$',alpha,c_l,alpha_ref_interp,cl_ref_interp,data.alpha_exp_cl,data.cl_exp,data.alpha_mod_cl,data.cl_mod,data.alpha_cycles,data.cl_cycles,authors,iti,source,[],model_name,plot_opt);
    [~,handles(2).fig] = alpha_plotter(0,handles(2).fig,'cm','$c_m$',alpha,c_m,alpha_ref_interp,cm_ref_interp,data.alpha_exp_cm,data.cm_exp,data.alpha_mod_cm,data.cm_mod,data.alpha_cycles,data.cm_cycles,authors,iti,source,[],model_name,plot_opt);
    [~,handles(3).fig] = alpha_plotter(0,handles(3).fig,'cd','$c_d$',alpha,c_d,alpha_ref_interp,cd_ref_interp,data.alpha_exp_cd,data.cd_exp,data.alpha_mod_cd,data.cd_mod,data.alpha_cycles,data.cd_cycles,authors,iti,source,[],model_name,plot_opt);
    [~,handles(4).fig] = alpha_plotter(0,handles(4).fig,'cn','$c_n$',alpha,c_n,alpha_ref_interp,cn_ref_interp,data.alpha_exp_cn,data.cn_exp,data.alpha_mod_cn,data.cn_mod,data.alpha_cycles,data.cn_cycles,authors,iti,source,[],model_name,plot_opt);
    [~,handles(5).fig] = alpha_plotter(0,handles(5).fig,'cc','$c_c$',alpha,c_c,alpha_ref_interp,cc_ref_interp,data.alpha_exp_cc,data.cc_exp,data.alpha_mod_cc,data.cc_mod,data.alpha_cycles,data.cc_cycles,authors,iti,source,[],model_name,plot_opt);
end

%% Coefficients x time
[~,handles(6).fig] = time_plotter(0,handles(6).fig,'cl x t','$c_l$',t,c_l,time_ref_interp,cl_ref_interp,data.time_exp_cl,data.clt_exp,data.time_mod_cl,data.clt_mod,data.time_cycles,data.cl_cycles,model,outputs,params,interp,authors,source,[],model_name,plot_opt);
[~,handles(7).fig] = time_plotter(0,handles(7).fig,'cm x t','$c_m$',t,c_m,time_ref_interp,cm_ref_interp,data.time_exp_cm,data.cmt_exp,data.time_mod_cm,data.cmt_mod,data.time_cycles,data.cm_cycles,model,outputs,params,interp,authors,source,[],model_name,plot_opt);
[~,handles(8).fig] = time_plotter(0,handles(8).fig,'cd x t','$c_d$',t,c_d,time_ref_interp,cd_ref_interp,data.time_exp_cd,data.cdt_exp,data.time_mod_cd,data.cdt_mod,data.time_cycles,data.cd_cycles,model,outputs,params,interp,authors,source,[],model_name,plot_opt);
[~,handles(9).fig] = time_plotter(0,handles(9).fig,'cn x t','$c_n$',t,c_n,time_ref_interp,cn_ref_interp,data.time_exp_cn,data.cnt_exp,data.time_mod_cn,data.cnt_mod,data.time_cycles,data.cn_cycles,model,outputs,params,interp,authors,source,[],model_name,plot_opt);
[~,handles(10).fig] = time_plotter(0,handles(10).fig,'cc x t','$c_c$',t,c_c,time_ref_interp,cc_ref_interp,data.time_exp_cc,data.cct_exp,data.time_mod_cc,data.cct_mod,data.time_cycles,data.cc_cycles,model,outputs,params,interp,authors,source,[],model_name,plot_opt);

%% Coefficients x plunge-induced AoA
if isfield(outputs,'alpha_plunge') && any(diff(plot_opt.xlim_alpha_plunge_vec))
    for j=11:15
        [handles(j),plot_opt] = init_plotter(struct,base_name,case_now,source,params,0);
    end
    alpha_plunge = outputs.alpha_plunge;
    [~,handles(11).fig] = alpha_plunge_plotter(0,handles(11).fig,'cl','$c_l$',alpha_plunge,c_l,data.alpha_exp_cl,data.cl_exp,data.alpha_mod_cl,data.cl_mod,authors,iti,[],model_name,plot_opt);
    [~,handles(12).fig] = alpha_plunge_plotter(0,handles(12).fig,'cm','$c_m$',alpha_plunge,c_m,data.alpha_exp_cm,data.cm_exp,data.alpha_mod_cm,data.cm_mod,authors,iti,[],model_name,plot_opt);
    [~,handles(13).fig] = alpha_plunge_plotter(0,handles(13).fig,'cd','$c_d$',alpha_plunge,c_d,data.alpha_exp_cd,data.cd_exp,data.alpha_mod_cd,data.cd_mod,authors,iti,[],model_name,plot_opt);
    [~,handles(14).fig] = alpha_plunge_plotter(0,handles(14).fig,'cn','$c_n$',alpha_plunge,c_n,data.alpha_exp_cn,data.cn_exp,data.alpha_mod_cn,data.cn_mod,authors,iti,[],model_name,plot_opt);
    [~,handles(15).fig] = alpha_plunge_plotter(0,handles(15).fig,'cc','$c_c$',alpha_plunge,c_c,data.alpha_exp_cc,data.cc_exp,data.alpha_mod_cc,data.cc_mod,authors,iti,[],model_name,plot_opt);
end

% %% Coefficients x flap angle
% if contains(model,["flap","gen"])
%     for j=11:16
%         [handles(j),plot_opt] = init_plotter(struct,base_name,case_now,source,params,0);
%     end
%     delta = outputs.delta;
%     c_h = outputs.c_h;
%     [~,handles(11).fig] = delta_plotter(0,handles(11).fig,'cl','$c_l$',delta,c_l,data.delta_exp_cl,data.cl_exp,data.delta_mod_cl,data.cl_mod,authors,iti,[],model_name,plot_opt);
%     [~,handles(12).fig] = delta_plotter(0,handles(12).fig,'cm','$c_m$',delta,c_m,data.delta_exp_cm,data.cm_exp,data.delta_mod_cm,data.cm_mod,authors,iti,[],model_name,plot_opt);
%     [~,handles(13).fig] = delta_plotter(0,handles(13).fig,'cd','$c_d$',delta,c_d,data.delta_exp_cd,data.cd_exp,data.delta_mod_cd,data.cd_mod,authors,iti,[],model_name,plot_opt);
%     [~,handles(14).fig] = delta_plotter(0,handles(14).fig,'cn','$c_n$',delta,c_n,data.delta_exp_cn,data.cn_exp,data.delta_mod_cn,data.cn_mod,authors,iti,[],model_name,plot_opt);
%     [~,handles(15).fig] = delta_plotter(0,handles(15).fig,'cc','$c_c$',delta,c_c,data.delta_exp_cc,data.cc_exp,data.delta_mod_cc,data.cc_mod,authors,iti,[],model_name,plot_opt);
%     [~,handles(16).fig] = delta_plotter(0,handles(16).fig,'ch','$c_h$',delta,c_h,data.delta_exp_ch,data.ch_exp,data.delta_mod_ch,data.ch_mod,authors,iti,[],model_name,plot_opt);
% end

%% Save to filepaths
saveas(handles(1).fig,filepath_cl_a);
saveas(handles(2).fig,filepath_cm_a);
saveas(handles(3).fig,filepath_cd_a);
saveas(handles(4).fig,filepath_cn_a);
saveas(handles(5).fig,filepath_cc_a);
saveas(handles(6).fig,filepath_cl_t);
saveas(handles(7).fig,filepath_cm_t);
saveas(handles(8).fig,filepath_cd_t);
saveas(handles(9).fig,filepath_cn_t);
saveas(handles(10).fig,filepath_cc_t);

end