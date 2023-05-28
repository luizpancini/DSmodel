function xdot = BL_state_space(x,xdot,w_QS,alpha_QS,w_tqc,q,R,A,fprime_n,fprime_m,fprime_c,Tf_n,Tf_m,Tf_c,Ta,Ta_theta)

% Potential flow states 
xdot(1) = A(1)*x(1)+w_tqc; 
xdot(2) = A(2)*x(2)+w_tqc; 
xdot(3) = A(3)*x(3)+alpha_QS;
xdot(4) = A(4)*x(4)+q;
xdot(5) = A(5)*x(5)+alpha_QS; 
xdot(6) = A(6)*x(6)+alpha_QS; 
xdot(7) = A(7)*x(7)+q; 
xdot(8) = A(8)*x(8)+q; 

% Nonlinear states
xdot(9)  = (w_QS-x(9))/Ta;              % Delayed quasi-steady pitch-plunge-induced downwash at pitch/plunge axis 
xdot(10) = (fprime_n-x(10))/Tf_n;       % Delayed separation point - normal force
xdot(11) = (fprime_m-x(11))/Tf_m;       % Delayed separation point - pitch moment
xdot(12) = (fprime_c-x(12))/Tf_c;       % Delayed separation point - chordwise force
xdot(13) = (R-x(13))/(3*Ta);            % Delayed capped reduced pitch rate
xdot(14) = (R-x(14))/Ta_theta;          % Delayed capped reduced pitch rate for the delayed AoA - increased lag under stall suppression conditions at relatively high pitch rates

end