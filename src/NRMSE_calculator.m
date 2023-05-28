function [cn_NRMSE,cm_NRMSE,cc_NRMSE] = NRMSE_calculator(c_n,c_m,c_c,cn_ref_interp,cm_ref_interp,cc_interp)

% Root mean square error normalized by the range
cn_NRMSE = sqrt(mean((c_n(:)-cn_ref_interp(:)).^2))/(max(cn_ref_interp)-min(cn_ref_interp));
cm_NRMSE = sqrt(mean((c_m(:)-cm_ref_interp(:)).^2))/(max(cm_ref_interp)-min(cm_ref_interp));
cc_NRMSE = sqrt(mean((c_c(:)-cc_interp(:)).^2))/(max(cc_interp)-min(cc_interp));

end