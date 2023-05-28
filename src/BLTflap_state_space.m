function xdot = BLTflap_state_space(x,xdot,alpha,alpha_tqc,delta_qs,R,A,Tf_n,Tf_m,Tf_c,Ta,fprime_n,fprime_m,fprime_c,Ta_theta)

% % Potential flow states from thin-airfoil theory
% xdot(1) = x(2); 
% xdot(2) = B*x(1:2)+U/b*alpha+alphadot;  

% Potential flow states: airfoil-induced and flap-induced
xdot(1)  = A(1)*x(1)+alpha_tqc; 
xdot(2)  = A(2)*x(2)+alpha_tqc; 
xdot(9)  = A(1)*x(9)+delta_qs; 
xdot(10) = A(2)*x(10)+delta_qs;

% Nonlinear states
xdot(3) = (alpha-x(3))/Ta;            % Delayed airfoil-induced AoA 
xdot(4) = (fprime_n-x(4))/Tf_n;       % Delayed separation point - normal force
xdot(5) = (fprime_m-x(5))/Tf_m;       % Delayed separation point - pitch moment
xdot(6) = (fprime_c-x(6))/Tf_c;       % Delayed separation point - chordwise force
xdot(7) = (R-x(7))/(Ta*3);            % Delayed reduced pitch rate
xdot(8) = (R-x(8))/(Ta_theta);        % Delayed reduced pitch rate for the delayed AoA - increased lag under stall suppression conditions at relatively high pitch rates

end