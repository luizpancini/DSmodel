function xdot = BLThargen_state_space(x,xdot,b,U,beta,w_QS,cnawdot_tqc,wdot_tqc_f,R,wg,Tf_n,Tf_m,Tf_c,Ta,fprime_n,fprime_m,fprime_c,Ta_theta,A1,A2,b1,b2,A1W,A2W,b1W,b2W,bG)

%% Nonlinear states
xdot(1) = (w_QS-x(1))/Ta;             % Delayed pitch-plunge-induced quasi-steady downwash 
xdot(2) = (fprime_n-x(2))/Tf_n;       % Delayed separation point - normal force
xdot(3) = (fprime_m-x(3))/Tf_m;       % Delayed separation point - pitch moment
xdot(4) = (fprime_c-x(4))/Tf_c;       % Delayed separation point - chordwise force
xdot(5) = (R-x(5))/(Ta*3);            % Delayed capped reduced pitch rate
xdot(6) = (R-x(6))/(Ta_theta);        % Decreased-lag delayed capped reduced pitch rate 
xdot(7) = (wg-x(7))/(5*Ta);           % Delayed gust-induced downwash

%% Potential flow states
% Pitch-plunge-flap-induced
xdot(8)  = -U/b*b1W*beta^2*x(8)+A1W*cnawdot_tqc;        % Circulatory pitch-plunge-pitch-rate-induced 3/4-chord downwash by c_n_alpha product state component
xdot(9)  = -U/b*b2W*beta^2*x(9)+A2W*cnawdot_tqc;        % Circulatory pitch-plunge-pitch-rate-induced 3/4-chord downwash by c_n_alpha product state component
xdot(10) = -U/b*b1W*beta^2*x(10)+A1W*wdot_tqc_f;        % Circulatory flap-induced 3/4-chord downwash state component
xdot(11) = -U/b*b2W*beta^2*x(11)+A2W*wdot_tqc_f;        % Circulatory flap-induced 3/4-chord downwash state component
% Gust-induced
xdot(12:end) = -U/b*beta^2*bG.*x(12:end)+wg;            % Gust-induced velocity state components

end