function outputs = BLS_output_vars(t,x,y,xdot)

outputs.x = x;
outputs.t = t;
outputs.alpha = y(1,:);
outputs.c_n = y(2,:);
outputs.c_m = y(3,:);
outputs.c_c = y(4,:);
outputs.c_l = outputs.c_n.*cos(outputs.alpha)+outputs.c_c.*sin(outputs.alpha);
outputs.c_d = outputs.c_n.*sin(outputs.alpha)-outputs.c_c.*cos(outputs.alpha);
outputs.f = y(5,:);
outputs.fprime = y(6,:);
outputs.f2prime = x(8,:);
outputs.c_nf = y(7,:);
outputs.c_nI = y(8,:);
outputs.c_nv = y(9,:);
outputs.c_mf = y(10,:);
outputs.c_mv = y(11,:);
outputs.c_mI = y(12,:);
outputs.tau_v = y(13,:);
outputs.dCP = y(14,:);
outputs.so_lim = y(15,:);
outputs.qR = y(16,:); 
outputs.alpha1 = y(17,:);
outputs.R = y(18,:);
outputs.Tf = y(19,:);
outputs.c_nC = y(20,:);
outputs.theta = x(7,:)./outputs.so_lim;

end