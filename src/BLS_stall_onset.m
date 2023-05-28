function [stall_state,alpha_cr,theta] = BLS_stall_onset(stall_state,R,alpha_ss,alpha_ds0)

% Critical angle for dynamic stall
alpha_cr = alpha_ss + (alpha_ds0-alpha_ss)*R;

% Dynamic stall onset ratio
theta = stall_state/alpha_cr;

end