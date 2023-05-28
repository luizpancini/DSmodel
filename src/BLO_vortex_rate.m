function c_vdot = BLO_vortex_rate(alpha,theta,tau_v,TvL,fprime,f2prime,Tf,c_nC,c_nCdot)

% Initialize 
c_vdot = 0;

% Check if vortex is over the chord
if theta>=1 && tau_v<=TvL 
    % Kirchhoff/Helmholtz factor and its rate
    K_f = ((1+sqrt(f2prime))/2)^2;
    K_fdot = 1/4*(1+1/sqrt(f2prime))*(fprime-f2prime)/Tf;
    % Vorticity coefficient rate
    c_vdot = c_nCdot*(1-K_f) - c_nC*K_fdot;
    % Check motion 
    if alpha*c_vdot < 0 % This is a suggestion by Bjorck, so the vortex feedback cannot be negative
        c_vdot = 0;
    end
end

end