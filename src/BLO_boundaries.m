function boundary = BLO_boundaries(boundary,t_ip1,t_i,tv0,alpha_i,alpha1_i,f2prime_i,q_i,so_i,alpha1_ip1,alpha_ip1,f2prime_ip1,q_ip1,so_ip1,alpha1_0,c_n1,c_n_alpha,fb,TvL)

% Discontinuity boundaries
boundary(1) = (abs(alpha_ip1)-alpha1_0)*(abs(alpha_i)-alpha1_0);        % f: alpha = alpha1_0
boundary(2) = (so_ip1-c_n1)*(so_i-c_n1);                                % so_state = so_lim
boundary(3) = (so_ip1/c_n_alpha-alpha1_i)*(so_i/c_n_alpha-alpha1_i);    % f': x9/c_n_alpha = alpha_1
boundary(4) = (t_ip1-tv0-TvL)*(t_i-tv0-TvL);                            % t = tv0+TvL
boundary(5) = (t_ip1-tv0-2*TvL)*(t_i-tv0-2*TvL);                        % t = tv0+2*TvL
boundary(6) = (alpha_ip1*q_ip1)*(alpha_i*q_i);                          % alpha*q = 0
boundary(8) = (f2prime_ip1-fb)*(f2prime_i-fb);                          % f2prime = fb

end