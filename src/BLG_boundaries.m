function boundary = BLG_boundaries(boundary,t_ip1,t_i,tv0,alpha_i,alpha1_i,alpha2_i,q_i,so_i,alpha1_ip1,alpha_ip1,alpha2_ip1,q_ip1,so_ip1,alpha1_0,alpha2_0,c_n1,c_n_alpha,TvL)

% Discontinuity boundaries
boundary(1) = (abs(alpha_ip1)-alpha1_0)*(abs(alpha_i)-alpha1_0);        % f: alpha = alpha1_0
boundary(2) = (abs(alpha_ip1)-alpha2_0)*(abs(alpha_i)-alpha2_0);        % f: alpha = alpha2_0
boundary(3) = (so_ip1/c_n_alpha-alpha1_i)*(so_i/c_n_alpha-alpha1_i);    % f': x9/c_n_alpha = alpha_1
boundary(4) = (so_ip1/c_n_alpha-alpha2_i)*(so_i/c_n_alpha-alpha2_i);    % f': x9/c_n_alpha = alpha_2
boundary(5) = (so_ip1-c_n1)*(so_i-c_n1);                                % so_state = so_lim
boundary(6) = (t_ip1-tv0-TvL)*(t_i-tv0-TvL);                            % t = tv0+TvL
boundary(7) = (t_ip1-tv0-2*TvL)*(t_i-tv0-2*TvL);                        % t = tv0+2*TvL
boundary(8) = (alpha_ip1*q_ip1)*(alpha_i*q_i);                          % alpha*q = 0

end