function outputs = BLO_output_vars(t,x,y,xdot)

outputs.x = x;
outputs.t = t;

outputs.alpha = y(1,:);
outputs.alpha1 = y(2,:);
outputs.alpha_E = y(3,:);
outputs.theta = y(4,:);
outputs.tau_v = y(5,:);
outputs.c_c = y(6,:);
outputs.c_m = y(7,:);
outputs.c_mC = y(8,:);
outputs.c_mf = y(9,:);
outputs.c_mI = y(10,:);
outputs.c_mv = y(11,:);
outputs.c_n = y(12,:);
outputs.c_nC = y(13,:);
outputs.c_nf = y(14,:);
outputs.c_nI = y(15,:);
outputs.c_nP = y(16,:);
outputs.c_nv = y(17,:);
outputs.dalpha1 = y(18,:);
outputs.dCP = y(19,:);
outputs.f = y(20,:);
outputs.fprime = y(21,:);
outputs.q = y(22,:);
outputs.Tf = y(23,:);
outputs.Tv = y(24,:);

outputs.f2prime = x(10,:);

outputs.c_l = outputs.c_n.*cos(outputs.alpha)+outputs.c_c.*sin(outputs.alpha);
outputs.c_d = outputs.c_n.*sin(outputs.alpha)-outputs.c_c.*cos(outputs.alpha);

end