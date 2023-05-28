function xdot = BLgust_state_space(x,xdot,alpha,alpha_tqc,q,R,A,Tf_n,Tf_m,Tf_c,Ta,fprime_n,fprime_m,fprime_c,Ta_theta,Ag,Bg,gust_states,alpha_g,U,b,beta,g1,g2)

% Potential flow states 
xdot(1) = A(1)*x(1)+alpha_tqc; 
xdot(2) = A(2)*x(2)+alpha_tqc; 
xdot(3) = A(3)*x(3)+alpha;
xdot(4) = A(4)*x(4)+q;
xdot(5) = A(5)*x(5)+alpha; 
xdot(6) = A(6)*x(6)+alpha; 
xdot(7) = A(7)*x(7)+q; 
xdot(8) = A(8)*x(8)+q; 

% Nonlinear states
xdot(9)  = (alpha-x(9))/Ta;                         % Delayed pitch-induced AoA 
xdot(10) = (fprime_n-x(10))/Tf_n;                   % Delayed separation point - normal force
xdot(11) = (fprime_m-x(11))/Tf_m;                   % Delayed separation point - pitch moment
xdot(12) = (fprime_c-x(12))/Tf_c;                   % Delayed separation point - chordwise force
xdot(13) = (R-x(13))/(Ta*3);                        % Delayed capped reduced pitch rate
xdot(14) = (R-x(14))/(Ta_theta);                    % Delayed capped reduced pitch rate for the delayed AoA - increased lag under stall suppression conditions at relatively high pitch rates
xdot(end)= (alpha_g-x(end))/(5/2*Ta);               % Delayed gust-induced AoA

% Gust states
xdot(15:end-3) = Ag*gust_states+Bg*alpha_g;         % Gust states
xdot(end-2)    = -U/b*beta^2*g1*x(end-2)+alpha_g;   % Gust-convection state
xdot(end-1)    = -U/b*beta^2*g2*x(end-1)+alpha_g;   % Gust-convection state

end