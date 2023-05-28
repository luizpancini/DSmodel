

a_0 = OUTPUTS.params.a_0;
a_1 = OUTPUTS.params.a_1;
b = OUTPUTS.params.b;
U = OUTPUTS.params.U;
k = OUTPUTS.params.k;
r0 = OUTPUTS.params.r0;
t_cycle = OUTPUTS.params.t_cycle;

t = 0:1e-4:2*t_cycle;
alpha = a_0+a_1*sin(k*U/b*t); sa = sin(alpha); ca = cos(alpha);
alphadot = a_1*k*U/b*cos(k*U/b*t);
w = U*sa; u_c = U*ca;
wdot = -U*ca.*alphadot; udot_c = U*sa.*alphadot;
figure; plot(t,w,t,u_c-U); grid

u_ip = sqrt(u_c.^2+w.^2);
alphadot_QS = (u_c.*wdot-w.*udot_c)./u_ip.^2;
r = alphadot_QS*b./u_ip;
qR = abs(r)/r0;
R = min(1, qR);

figure; plot(alpha*180/pi,qR/b); grid