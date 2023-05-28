function c_c = BLO_cc_coeff(alpha,so_state,f2prime,theta,alpha_E,alpha_0L,eta,c_d0,c_n1,c_n_alpha,Df,E0)
 
if theta < 1
    c_c = -c_d0*cos(alpha)+eta*c_n_alpha*(alpha_E-alpha_0L)*sin(alpha_E-alpha_0L)*(sqrt(f2prime)-E0); 
else
    phi = Df*(abs(so_state)-c_n1); 
    c_c = -c_d0*cos(alpha)+eta*c_n_alpha*(alpha_E-alpha_0L)*sin(alpha_E-alpha_0L)*(sqrt(f2prime)*f2prime^phi-E0); % See Leishman and Beddoes (1986) - A Generalised Model... - page 11
end

end

