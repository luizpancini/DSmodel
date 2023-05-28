function [c_nC,c_nP,c_nCdot] = BLO_cn_vars(x,alpha,q,U,b,M,beta,c_n_alpha,A1,A2,b1,b2,K_a,K_q,T_I)

% Circulatory normal coefficient
c_nC = c_n_alpha*beta^2*U/b*(A1*b1*x(1)+A2*b2*x(2));

% Circulatory normal coefficient rate
c_nCdot = c_n_alpha*U/b*beta^2*(A1*b1*(-U/b*beta^2*b1*x(1)+alpha+q/2)+A2*b2*(-U/b*beta^2*b2*x(2)+alpha+q/2));

% Impulsive normal coefficient
c_nI = - 4/(M*K_a*T_I)*x(3)-1/(M*K_q*T_I)*x(4)+4/M*alpha+1/M*q;

% Total potential flow normal coefficient
c_nP = c_nC + c_nI;

end