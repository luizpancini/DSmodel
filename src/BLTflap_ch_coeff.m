function [c_h,c_hC,c_hI] = BLTflap_ch_coeff(alphadot,alphaddot,delta,deltadot,deltaddot,T1,T3,T4,T5,T9,T10,T11,T12,T13,eps_h,alpha_0L,U,b,alpha_Et)

% Circulatory unsteady
c_hC = -eps_h*T12/2*(alpha_Et-alpha_0L); 

% Quasi-steady + inertial: airfoil-induced, flap-induced and total
c_hIa = -1/(2*U^2)*(U*b*(-2*T9-T1-T4)*alphadot+2*T13*b^2*alphaddot);
c_hIf = -1/(2*U^2)*(U^2/pi*(T5-T4*T10)*delta-U*b/2/pi*T4*T11*deltadot-T3/pi*b^2*deltaddot);
c_hI = 1*(c_hIa+c_hIf);

% Total: circulatory + quasi-steady + inertial 
c_h = c_hC+c_hI;

end
