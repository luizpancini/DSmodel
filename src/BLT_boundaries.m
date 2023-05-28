function boundary = BLT_boundaries(t_ip1,t_i,x10_ip1,x10_i,x12_ip1,x12_i,tv0,so_i,so_lim_i,alpha1_i,alpha_i,alpha_dot_i,q_i,so_ip1,so_lim_ip1,alpha1_ip1,alpha_ip1,alpha_dot_ip1,q_ip1,alpha1_0,TvL,q0)

% Discontinuity boundaries
boundary(1)=(abs(alpha_ip1)-alpha1_0)*(abs(alpha_i)-alpha1_0);        % f: alpha and alpha1_0
boundary(2)=(so_ip1-so_lim_ip1)*(so_i-so_lim_i);                      % so_state and so_lim
boundary(3)=(so_ip1-alpha1_i)*(so_i-alpha1_i);                        % f': x9 and alpha_1
boundary(4)=(t_ip1-tv0-TvL)*(t_i-tv0-TvL);
boundary(5)=(t_ip1-tv0-2*TvL)*(t_i-tv0-2*TvL);
boundary(6)=(alpha_ip1*alpha_dot_ip1)*(alpha_i*alpha_dot_i);
boundary(7)=(x10_ip1-x12_ip1)*(x10_i-x12_i);
boundary(8)=(x10_ip1-0.7)*(x10_i-0.7);
boundary(9) = (abs(q_ip1)-q0)*(abs(q_i)-q0);
boundary(10) = q_ip1*q_i;

end