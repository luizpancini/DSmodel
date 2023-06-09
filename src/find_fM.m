function y = find_fM(f2prime_m,c_mf,c_nf,K0,K1,K2,K3,kappa,R,RD,S,P,theta,q)

if theta*q > 0
    d = -R*abs(theta)*(1-f2prime_m);
else
    d = 0; 
    K2 = K2*(1+2*R*RD*S*(1-P)); 
end
y = (K0+K1*(1-f2prime_m)+K2*sin(pi*f2prime_m^kappa)+K3*d)*c_nf-c_mf;

end