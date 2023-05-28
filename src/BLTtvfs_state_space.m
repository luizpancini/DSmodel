function xdot = BLTtvfs_state_space(x,xdot,alpha,R,Tf_n,Tf_m,Tf_c,Ta,fprime_n,fprime_m,fprime_c,Ta_theta,U,b,A1,A2,b1,b2,w_tqc_dot)

% Nonlinear states
xdot(1) = (alpha-x(1))/Ta;            % Delayed AoA 
xdot(2) = (fprime_n-x(2))/Tf_n;       % Delayed separation point - normal force
xdot(3) = (fprime_m-x(3))/Tf_m;       % Delayed separation point - pitch moment
xdot(4) = (fprime_c-x(4))/Tf_c;       % Delayed separation point - chordwise force
xdot(5) = (R-x(5))/(Ta*3);            % Delayed reduced pitch rate
xdot(6) = (R-x(6))/(Ta_theta);        % Delayed reduced pitch rate for the delayed AoA - increased lag under stall suppression conditions at relatively high pitch rates

% Potential flow states
xdot(7) = -U/b*b1*x(7)+A1*w_tqc_dot; 
xdot(8) = -U/b*b2*x(8)+A2*w_tqc_dot;

end