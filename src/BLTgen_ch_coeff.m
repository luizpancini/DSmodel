function [c_h,c_hC,c_hI] = BLTgen_ch_coeff(U,b,ah,w_Ep,w_Ef,alphadot,alphaddot,hddot,delta,deltadot,deltaddot,T1,T3,T4,T5,T9,T10,T11,T12,T13,alpha_0L,eps_fh)

% Circulatory unsteady 
c_hC = -T12/(2*U)*(w_Ep-U*sin(alpha_0L)+eps_fh*w_Ef);

% Quasi-steady
c_hqs = -1/(2*U^2)*((-2*T9-T1+T4*(ah-1/2))*U*b*alphadot+U^2/pi*(T5-T4*T10)*delta-U*b/2/pi*T4*T11*deltadot);

% Inertial: pitch-plunge-induced, flap-induced and total
c_hIa = -1/(2*U^2)*(2*T13*b^2*alphaddot-T1*b*hddot);
c_hIf = -1/(2*U^2)*(-T3/pi*b^2*deltaddot);
c_hI = c_hIa+c_hIf;

% Total: circulatory + quasi-steady + inertial 
c_h = c_hC+c_hqs+c_hI;

end
