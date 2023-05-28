function outputs = BLgust_output_vars(t,x,y,xdot)

outputs.x = x;
outputs.t = t;

outputs.alpha = y(1,:);
outputs.c_n = y(2,:);
outputs.c_m = y(3,:);
outputs.c_c = y(4,:);
outputs.c_l = outputs.c_n.*cos(outputs.alpha)+outputs.c_c.*sin(outputs.alpha);
outputs.c_d = outputs.c_n.*sin(outputs.alpha)-outputs.c_c.*cos(outputs.alpha);
outputs.c_nf = y(5,:);
outputs.c_nI = y(6,:);
outputs.c_mf = y(7,:);
outputs.c_mI = y(8,:);
outputs.c_nv = y(9,:);
outputs.c_mv = y(10,:);
outputs.c_cv = y(11,:);
outputs.c_lv = outputs.c_nv.*cos(outputs.alpha)+outputs.c_cv.*sin(outputs.alpha);
outputs.f = y(12,:);
outputs.tau_v = y(13,:);
outputs.dCP = y(14,:);
outputs.so_lim = y(15,:);
outputs.qR = y(16,:); 
outputs.alpha1_n = y(17,:);
outputs.K_f = y(18,:);
outputs.c_nC = y(19,:);
outputs.Tf_n = y(20,:);
outputs.dalpha1_n = y(21,:);
outputs.fprime_n = y(22,:);
outputs.fprime_m = y(23,:);
outputs.fprime_c = y(24,:);
outputs.q = y(25,:);
outputs.dalpha1_m = y(26,:);
outputs.dalpha1_c = y(27,:);
outputs.R_dot = y(28,:);
outputs.Tf_m = y(29,:);
outputs.Tf_c = y(30,:);
outputs.theta_min = y(31,:);
outputs.theta_max = y(32,:);
outputs.P = y(33,:);
outputs.S = y(34,:);
outputs.alpha_E = y(35,:);
outputs.alpha_g = y(36,:);
outputs.R = outputs.qR; dummy = (outputs.R>1); outputs.R(dummy) = 1;

outputs.so_state = x(9,:)+x(end,:);
outputs.theta = outputs.so_state./outputs.so_lim;
outputs.f2prime_n = x(10,:);
outputs.f2prime_m = x(11,:);
outputs.f2prime_c = x(12,:);
outputs.RD = x(13,:);
outputs.RD_theta = x(14,:);

end