function [alpha1,alpha2,dalpha1,dalpha2] =  BLG_alpha_brk(alpha,q,alpha1_0,alpha2_0,delta_alpha1,delta_alpha2,f2prime)

% Breakpoint of separation angles offsets
dalpha1 = 0; dalpha2 = 0;
if alpha*q<0
    dalpha1 = -(1-f2prime)^(1/4)*delta_alpha1;
    dalpha2 = -(1-f2prime)^(1/4)*delta_alpha2;
end

% Total unsteady breakpoint of separation angles
alpha1 = alpha1_0+dalpha1; 
alpha2 = alpha2_0+dalpha2; 

end