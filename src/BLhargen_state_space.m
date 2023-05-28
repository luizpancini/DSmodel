function xdot = BLhargen_state_space(x,xdot,b,ah,U,M,Mdot,beta,w_QS,alpha_QS,alphadot_QS,alphaddot,cnawdot_tqc,wdot_tqc_f,R,wg,Ucm_QS_dot,Tf_n,Tf_m,Tf_c,Ta,fprime_n,fprime_m,fprime_c,Ta_theta,T_a,T_q,T_M,T_am,T_qm,T_Mm,A1,A2,A3,A4,b1,b2,b3,b4,b5,A1W,A2W,b1W,b2W,bG)

%% Nonlinear states
xdot(1) = (w_QS-x(1))/Ta;             % Delayed quasi-steady pitch-plunge-induced downwash 
xdot(2) = (fprime_n-x(2))/Tf_n;       % Delayed separation point - normal force
xdot(3) = (fprime_m-x(3))/Tf_m;       % Delayed separation point - pitch moment
xdot(4) = (fprime_c-x(4))/Tf_c;       % Delayed separation point - chordwise force
xdot(5) = (R-x(5))/(Ta*3);            % Delayed capped reduced pitch rate
xdot(6) = (R-x(6))/(Ta_theta);        % Decreased-lag delayed capped reduced pitch rate 
xdot(7) = (wg-x(7))/(5*Ta);           % Delayed gust-induced downwash

%% Potential flow states
% Airfoil motion-induced
xdot(8)  = -U/b*b1*beta^2*x(8)+A1*cnawdot_tqc;              % Circulatory effective pitch-plunge-pitch-rate-induced 3/4-chord downwash by c_n_alpha product state component. Could also use A1*cnawdot_tqc*exp(-b1*ds/2) for increased accuracy (see Leishman Eq 8.97). The alternate state equation would be xdot(8) = -U/b*b1*beta^2*x(8)+w_tqc (considering constant c_n_alpha)
xdot(9)  = -U/b*b2*beta^2*x(9)+A2*cnawdot_tqc;              % Circulatory effective pitch-plunge-pitch-rate-induced 3/4-chord downwash by c_n_alpha product state component. The alternate state equation would be xdot(9) = -U/b*b2*beta^2*x(9)+w_tqc (considering constant c_n_alpha)
xdot(10) = -U/b*b1W*beta^2*x(10)+A1W*wdot_tqc_f;            % Circulatory effective flap-induced 3/4-chord downwash state component
xdot(11) = -U/b*b2W*beta^2*x(11)+A2W*wdot_tqc_f;            % Circulatory effective flap-induced 3/4-chord downwash state component
xdot(12) = -U/b/T_a*x(12)+4/M*alphadot_QS;                  % Non-circulatory pitch-plunge-induced c_n state component
xdot(13) = -U/b/T_q*x(13)-4*b*ah/(U*M)*alphaddot;           % Non-circulatory pitch-rate-induced c_n state component
xdot(14) = -U/b/T_M*x(14)+4*alpha_QS/M^2*Mdot;              % Non-circulatory airspeed-rate-induced c_n state component
xdot(15) = -U/b/(b3*T_am)*x(15)+2*ah*A3*1/M*alphadot_QS;    % Non-circulatory pitch-plunge-induced pitch axis c_m state component       
xdot(16) = -U/b/(b4*T_am)*x(16)+2*ah*A4*1/M*alphadot_QS;    % Non-circulatory pitch-plunge-induced pitch axis c_m state component 
xdot(17) = -U/b/(b3*T_Mm)*x(17)+2*ah*A3*alpha_QS/M^2*Mdot;  % Non-circulatory airspeed-rate-induced pitch axis c_m state component
xdot(18) = -U/b/(b4*T_Mm)*x(18)+2*ah*A4*alpha_QS/M^2*Mdot;  % Non-circulatory airspeed-rate-induced pitch axis c_m state component
xdot(19) = -U/b/T_qm*x(19)-(ah^2+1/3)*2*b/(U*M)*alphaddot;  % Non-circulatory pitch-rate-induced pitch axis c_m state component
xdot(20) = -U/b*beta^2*b5*x(20)+Ucm_QS_dot;                 % Circulatory (quasi-steady) pitch-rate-induced pitch axis c_m state component
% Gust-induced
xdot(21:end) = -U/b*beta^2*bG.*x(21:end)+wg;                % Gust-induced velocity state components

end