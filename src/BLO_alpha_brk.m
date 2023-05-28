function [alpha1,dalpha1] =  BLO_alpha_brk(alpha,q,alpha1_0,delta_alpha1,f2prime)

% Breakpoint of separation angle offset
dalpha1 = 0;
if alpha*q<0
    dalpha1 = -(1-f2prime)^(1/4)*delta_alpha1;
end

% Total unsteady breakpoint of separation angle
alpha1 = alpha1_0+dalpha1; 

end