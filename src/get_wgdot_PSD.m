function wgdot = get_wgdot_PSD(t,wg,dt)

% Initialize
wgdot = zeros(length(t),1);

%% Employ finite differences to calculate the derivatives
% First point - 2nd order forward difference
wgdot(1) = (-wg(3)+4*wg(2)-3*wg(1))/(2*dt);
% Second point - 1st order centered difference
wgdot(2) = (wg(3)-wg(1))/(2*dt);
% Last but one - 1st order centered difference
wgdot(end-1) = (wg(end)-wg(end-2))/(2*dt);
% Last point - 2nd order backward difference
wgdot(end) = (3*wg(end)-4*wg(end-1)+wg(end-2))/(2*dt);
% Middle points - 2nd order centered difference
for i=3:length(t)-2
    wgdot(i) = (-wg(i+2)+8*wg(i+1)-8*wg(i-1)+wg(i-2))/(12*dt);
end

end