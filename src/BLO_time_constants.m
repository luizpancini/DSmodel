function [Tf,Tv] = BLO_time_constants(alpha,q,theta,tau_v,TvL,Tf0,Tv0,f2prime,fb)

if theta >= 1 % Vortex shedding phase
    % Primary vortex
    if tau_v<=TvL && alpha*q>=0
        Tf = Tf0*3/4;   % Accelerate the rate of boundary layer dettachment during the vortex convection
        Tv = Tv0;       % Nominal rate of vortex accumulation
    elseif tau_v>TvL && tau_v<=2*TvL && alpha*q>=0
        Tf = Tf0/2;     % Accelerate more the rate of dettachment after the vortex reaches the trailing edge 
        Tv = Tv0/2;     % Increase the rate of decay of the vortex after it reaches the trailing edge 
    elseif tau_v>2*TvL && alpha*q>=0
        Tf = Tf0/2;     % Maintain the rate of boundary layer dettachment after the vortex is totally shed
        Tv = Tv0/2;     % Maintain a high rate of vortex lift decay after the vortex is totally shed
    elseif tau_v<=2*TvL && alpha*q<0
        Tf = Tf0/2;     % Maintain the rate of dettachment if the rate of change of AoA changes during the vortex shedding 
        Tv = Tv0/2;     % Increase the rate of decay of the vortex lift if the rate of change of AoA changes during the vortex shedding 
    elseif tau_v>2*TvL && alpha*q<0
        Tf = 4*Tf0;     % Delay the reattachment of the boundary layer after the vortex is totally shed and the AoA is decreasing
        Tv = Tv0/2;     % Maintain a high rate of vortex lift decay after the vortex is totally shed
    end
else % Reattachment phase
    Tv = Tv0/2;     % Maintain high rate of vortex lift decay after the vortex is totally shed
    Tf = 4*Tf0;     % Delay the reattachment of the boundary layer
    if f2prime>=fb && alpha*q>=0 % Dimitriadis' suggestion 
        Tf = Tf0;   % Set to nominal conditions if the rate of change of AoA is increasing and the flow is lightly separated
    elseif f2prime<fb && alpha*q>=0
        Tf = Tf0/2; % Accelerate boundary layer reattachment if the rate of change of AoA is already increasing and the flow is still massively separated
    end
end

end