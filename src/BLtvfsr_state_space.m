function xdot = BLtvfsr_state_space(x,xdot,alpha,R,Tf_n,Tf_m,Tf_c,Ta,fprime_n,fprime_m,fprime_c,Ta_theta)

% Nonlinear states
xdot(1) = (alpha-x(1))/Ta;            % Delayed AoA 
xdot(2) = (fprime_n-x(2))/Tf_n;       % Delayed separation point - normal force
xdot(3) = (fprime_m-x(3))/Tf_m;       % Delayed separation point - pitch moment
xdot(4) = (fprime_c-x(4))/Tf_c;       % Delayed separation point - chordwise force
xdot(5) = (R-x(5))/(Ta*3);            % Delayed capped reduced pitch rate
xdot(6) = (R-x(6))/(Ta_theta);        % Delayed capped reduced pitch rate for the delayed AoA - increased lag under stall suppression conditions at relatively high pitch rates

end