function xdot = BLtvfs_state_space(x,xdot,U,M,b,beta,Mdot,alpha,alphadot,alphaddot,cnaw_tqc_dot,Ucmqs_dot,R,Tf_n,Tf_m,Tf_c,Ta,fprime_n,fprime_m,fprime_c,Ta_theta,A1,A2,A3,A4,b1,b2,b3,b4,b5)

% Nonlinear states
xdot(1) = (alpha-x(1))/Ta;            % Delayed AoA 
xdot(2) = (fprime_n-x(2))/Tf_n;       % Delayed separation point - normal force
xdot(3) = (fprime_m-x(3))/Tf_m;       % Delayed separation point - pitch moment
xdot(4) = (fprime_c-x(4))/Tf_c;       % Delayed separation point - chordwise force
xdot(5) = (R-x(5))/(Ta*3);            % Delayed capped reduced pitch rate
xdot(6) = (R-x(6))/(Ta_theta);        % Delayed capped reduced pitch rate for the delayed AoA - increased lag under stall suppression conditions at relatively high pitch rates

% Potential flow states
[T_a,T_q,T_M,T_am,T_qm,T_Mm] = get_inertial_ind_params(A1,A2,A3,A4,b1,b2,b3,b4,b5,M,beta);
xdot(7)  = -U/b*b1*beta^2*x(7)+A1*cnaw_tqc_dot; 
xdot(8)  = -U/b*b2*beta^2*x(8)+A2*cnaw_tqc_dot;
xdot(9)  = -U/b/T_a*x(9)+4/M*alphadot;
xdot(10) = -U/b/T_q*x(10)+2*b/(U*M)*alphaddot;                 % Notice that 1/M*qdot ~= 2*b/(U*M)*alphaddot
xdot(11) = -U/b/T_M*x(11)+4*alpha/M^2*Mdot;
xdot(12) = -U/b/(b3*T_am)*x(12)-A3*1/M*alphadot;
xdot(13) = -U/b/(b4*T_am)*x(13)-A4*1/M*alphadot;
xdot(14) = -U/b/(b3*T_Mm)*x(14)-A3*alpha/M^2*Mdot;
xdot(15) = -U/b/(b4*T_Mm)*x(15)-A4*alpha/M^2*Mdot;
xdot(16) = -U/b/T_qm*x(16)-7*2*b/(12*U*M)*alphaddot;           % Notice that 7/(12*M)*qdot ~= 7*2*b/(12*U*M)*alphaddot
xdot(17) = -U/b*beta^2*b5*x(17)+Ucmqs_dot;

end