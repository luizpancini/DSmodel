function xdot = BLTgust_state_space(x,xdot,alpha,alpha_tqc,alpha_g,R,A,Ag,Bg,Tf_n,Tf_m,Tf_c,Ta,fprime_n,fprime_m,fprime_c,Ta_theta,gust_states)

% Potential flow states
xdot(1) = A(1)*x(1)+alpha_tqc; 
xdot(2) = A(2)*x(2)+alpha_tqc; 

% Nonlinear states
xdot(3) = (alpha-x(3))/Ta;                  % Delayed motion-induced AoA 
xdot(4) = (fprime_n-x(4))/Tf_n;             % Delayed separation point - normal force
xdot(5) = (fprime_m-x(5))/Tf_m;             % Delayed separation point - pitch moment
xdot(6) = (fprime_c-x(6))/Tf_c;             % Delayed separation point - chordwise force
xdot(7) = (R-x(7))/(Ta*3);                  % Delayed reduced pitch rate
xdot(8) = (R-x(8))/(Ta_theta);              % Delayed reduced pitch rate for the delayed AoA - increased lag under stall suppression conditions at relatively high pitch rates

% Gust states
xdot(9:end-1) = Ag*gust_states+Bg*alpha_g;  % Linear gust states
xdot(end)  = (alpha_g-x(end))/(10/3*Ta);    % Delayed gust-induced AoA

end