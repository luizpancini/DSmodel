function c_c = BLG_cc_coeff(alpha,c_nP_prime,f,f2prime,f2prime_tv0,theta,alpha_E,alpha_0L,eta,c_d0,c_n1,c_n_alpha,Df,E0,Ef)

f_cc = f2prime;

if theta < 1
    c_c = -c_d0*cos(alpha)+eta*c_n_alpha*(alpha_E-alpha_0L)*sin(alpha_E-alpha_0L)*(sqrt(f_cc)-E0); 
else
    phi = Df*(abs(c_nP_prime)-c_n1)+0*Ef*(f2prime - f2prime_tv0); 
    c_c = -c_d0*cos(alpha)+eta*c_n_alpha*(alpha_E-alpha_0L)*sin(alpha_E-alpha_0L)*(sqrt(f_cc)*f_cc^phi-E0); 
end

end

