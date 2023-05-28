function [cn_RMSE,cm_RMSE,cc_RMSE] = RMSE_calculator(c_n,c_m,c_c,cn_ref_interp,cm_ref_interp,cc_interp)

% Root mean square error 
cn_RMSE = sqrt(mean((c_n(:)-cn_ref_interp(:)).^2));
cm_RMSE = sqrt(mean((c_m(:)-cm_ref_interp(:)).^2));
cc_RMSE = sqrt(mean((c_c(:)-cc_interp(:)).^2));

end