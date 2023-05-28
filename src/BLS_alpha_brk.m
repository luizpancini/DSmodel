function alpha1 = BLS_alpha_brk(alpha1_0,alpha_ss,alpha_ds0,R,alpha,q)

% Breakpoint angle offsets
dalpha1 = (alpha_ds0-alpha_ss)*abs(R)*sign(alpha*q);

% Unsteady breakpoint of separation angle
alpha1 = alpha1_0+dalpha1; 

end

