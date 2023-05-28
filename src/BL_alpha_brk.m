function [alpha1_n,dalpha1_n,alpha1_m,dalpha1_m,alpha1_c,dalpha1_c,dalpha1_db] = BL_alpha_brk(qR,R,RD,RD_theta,del_RD_acc,RD_m,theta,theta_min,theta_max,R_max,S,P,T,T_s,T_ub,T_db,upstroke,downstroke_beginning,dalpha1_n_i,dalpha1_db,alpha1_0n,alpha1_0m,alpha1_0c,alpha_ds0,alpha_ss,delta_alpha_0,delta_alpha_1,d_cm,d_cc,z_cm,z_ccd,z_ccu)
    
%% Breakpoint angle offsets
sqrt_of_R = sqrt(R); 
if upstroke 
    f_R = R/R_max*sqrt_of_R;                                      % Goes to zero at the end of the upstorke, less effect on low pitch rates
    dalpha1_n = (alpha_ds0-alpha_ss)*RD;           
    dalpha1_m = dalpha1_n+d_cm*f_R;                                         
    dalpha1_c = dalpha1_n*(1-z_ccu*del_RD_acc)+d_cc*f_R;                                    
else 
    dalpha1_n = -S*(delta_alpha_0+delta_alpha_1*qR*(1+P));    % S to apply only if stalled, *(1+P) to increase effect on light stall
    dalpha1_m = dalpha1_n*(1-z_cm*sqrt_of_R*(1-P/2));                                  
    dalpha1_c = dalpha1_n*(1-z_ccd*sqrt_of_R*(1-P/4));      
end

%% Unsteady breakpoint of separation angles
alpha1_n = alpha1_0n+dalpha1_n; 
alpha1_m = alpha1_0m+dalpha1_m;
alpha1_c = alpha1_0c+dalpha1_c;

end