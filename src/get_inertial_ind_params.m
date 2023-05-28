function [T_a,T_q,T_M,T_am,T_qm,T_Mm] = get_inertial_ind_params(A1,A2,A3,A4,b1,b2,b3,b4,b5,M,beta)

k_a = 0.75; k_am = 1.0; k_q = k_a; k_qm = k_am;

T_a = 2*M*k_a/(1-M+pi*beta*M^2*(A1*b1+A2*b2));
T_q = 2*M*k_q/(1-M+2*pi*beta*M^2*(A1*b1+A2*b2));
T_M = 2*M*k_a/(1-M+pi/beta*M^2*(A1*b1+A2*b2));
T_am = 2*M*k_am*(A3*b4+A4*b3)/(b3*b4*(1-M));
T_qm = 2*M*k_qm*7/(15*(1-M)+3*pi*beta*M^2*b5);
T_Mm = T_am;

end