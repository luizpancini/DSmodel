function [tv0,tau_v,f2prime_tv0] = stall_time_BLG(tv0,so_state,so_lim,so_i,so_lim_i,t,t_i,f2prime_tv0,f2prime,f)

% Primary vortex
if abs(so_state)>=so_lim && abs(so_i)<so_lim_i && t>t_i
    % Swap signs for next calculations in case of vortex at negative AoA
    if so_state < 0
        so_lim = -so_lim;
        so_lim_i = -so_lim_i;
    end
    % Linear interpolation for time of vortex shedding
    m_so = (so_state-so_i)/(t-t_i);
    m_so_lim = (so_lim-so_lim_i)/(t-t_i);
    tv0 = t_i + (so_lim_i-so_i)/(m_so-m_so_lim);
    % Delayed separation point at the time of vortex shedding
    f2prime_tv0 = f2prime;
end
% Time since primary vortex shedding 
tau_v = max([0, t-tv0]); 

end