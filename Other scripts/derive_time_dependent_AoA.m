clc
clear 

syms U(t) alpha(t) hdot(t) b ah dh delta(t) T4 T10 T1 T8 T11

Udot = diff(U);
alphadot = diff(alpha);
alphaddot = diff(alphadot);
hddot = diff(hdot);
deltadot = diff(delta);

alpha_bar = atan((U*sin(alpha)+hdot*cos(alpha)-b*(1/2+ah)*alphadot)/(U*cos(alpha)-hdot*sin(alpha)));
alphadot_bar = simplify(diff(alpha_bar))
pretty(alphadot_bar)
alphaddot_bar = simplify(diff(alphadot_bar))
pretty(alphaddot_bar)

alpha_h = atan(hdot/U);
alphaddot_h = simplify(diff(alpha_h,2))
pretty(alphaddot_h)

Ucm_qs = -1/(2*U)*(pi*b*U*(1/2-ah)*alphadot);
dUcm_qs = diff(Ucm_qs)
pretty(dUcm_qs)