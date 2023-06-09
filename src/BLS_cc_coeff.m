function c_c = BLS_cc_coeff(alpha,so_state,so_lim,f2prime,theta,c_n_alpha,alpha_E,E0,Df,alpha_0L,c_d0,f_prime,eta)
 
% Set f_cc
f_cc = min([f_prime; f2prime]);

% if theta < 1
    c_c = -c_d0*cos(alpha)+eta*c_n_alpha*(alpha_E-alpha_0L)^2*(sqrt(f_cc)-E0);
% else
%     PHI = Df*c_n_alpha*(so_state-so_lim);
%     c_c = -c_d0*cos(alpha)+eta*c_n_alpha*(alpha_E-alpha_0L)^2*(sqrt(f_cc)*f_cc^PHI-E0); 
% end

end

