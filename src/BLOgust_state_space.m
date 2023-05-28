function xdot = BLOgust_state_space(x,xdot,alpha,alpha_tqc,q,A,B,f,fprime,Tp,Tf,Tv,Tf0,c_nP,c_vdot,gust_states_rates)

% Potential flow states 
xdot(1) = A(1)*x(1)+alpha_tqc; 
xdot(2) = A(2)*x(2)+alpha_tqc; 
xdot(3) = A(3)*x(3)+B(3,1)*alpha+B(3,2)*q;
xdot(4) = A(4)*x(4)+B(4,1)*alpha+B(4,2)*q;
xdot(5) = A(5)*x(5)+B(5,1)*alpha+B(5,2)*q; 
xdot(6) = A(6)*x(6)+B(6,1)*alpha+B(6,2)*q; 
xdot(7) = A(7)*x(7)+B(7,1)*alpha+B(7,2)*q; 
xdot(8) = A(8)*x(8)+B(8,1)*alpha+B(8,2)*q; 

% Nonlinear states
xdot(9) = (c_nP-x(9))/Tp;           % Delayed potential flow normal coefficient 
xdot(10) = (fprime-x(10))/Tf;       % Delayed separation point 
xdot(11) = c_vdot-x(11)/Tv;         % Delayed rate of vortex accumulation
xdot(12) = (f-x(12))/(Tf0/2);       % Delayed separation point for pitching moment coefficient

% Gust states
xdot(13:end) = gust_states_rates;

end