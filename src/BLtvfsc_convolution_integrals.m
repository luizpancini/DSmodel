function [int_conv_Ia,int_conv_Iq,int_conv_IM,int_conv_Iam,int_conv_Iqm,int_conv_IMm,cnaw_tqc,w_tqc_eff,Ucm_qs_eff,s] = BLtvfsc_convolution_integrals(t,a_inf,b,U,M,beta,alpha,alphadot,q,w_tqc,Ucm_qs,c_n_alpha,A1,A2,A3,A4,b1,b2,b3,b4,b5,tvec,cnaw_tqc_vec,alpha_vec,alphadot_vec,M_vec,Ucm_qs_vec,s_vec,q_vec)

% Current inertial indicial parameters
[T_a,T_q,T_M,T_am,T_qm,T_Mm] = get_inertial_ind_params(A1,A2,A3,A4,b1,b2,b3,b4,b5,M,beta);

% Current circulatory variables 
cnaw_tqc = c_n_alpha*w_tqc;

% Arrays and their differences
cnaw_tqc_vec = [cnaw_tqc_vec; cnaw_tqc]; dcnaw_tqc = diff(cnaw_tqc_vec);
alpha_vec = [alpha_vec; alpha]; dalpha = diff(alpha_vec);
alphadot_vec = [alphadot_vec; alphadot]; dalphadot = diff(alphadot_vec);
M_vec = [M_vec; M]; dM = diff(M_vec);
Ucm_qs_vec = [Ucm_qs_vec; Ucm_qs]; dUcmqs = diff(Ucm_qs_vec);
q_vec = [q_vec; q]; dq = diff(q_vec);

% Semi-chords travelled: increment, array and current value
ds = a_inf*(M+M_vec(end-1))/(2*b)*(tvec(end)-tvec(end-1));
s_vec = [s_vec; s_vec(end)+ds];
s = s_vec(end);

% Range of indices for convolution
N = length(dM);
range = 1:N;

% Indicial functions
ind_fun_C = A1*exp(-b1*beta^2*(s-s_vec(range)))+A2*exp(-b2*beta^2*(s-s_vec(range)));
ind_fun_Cm = exp(-b5*beta^2*(s-s_vec(range)));
ind_fun_Ia = exp(-1/T_a*(s-s_vec(range)));
ind_fun_Iq = exp(-1/T_q*(s-s_vec(range)));
ind_fun_IM = exp(-1/T_M*(s-s_vec(range)));
ind_fun_Iam = A3*exp(-1/b3/T_am*(s-s_vec(range)))+A4*exp(-1/b4/T_am*(s-s_vec(range)));
ind_fun_Iqm = exp(-1/T_qm*(s-s_vec(range)));
ind_fun_IMm = A3*exp(-1/b3/T_Mm*(s-s_vec(range)))+A4*exp(-1/b4/T_Mm*(s-s_vec(range)));

% Convolution integrals
int_conv_C =  sum(dcnaw_tqc.*ind_fun_C);
int_conv_Cm = sum(dUcmqs.*ind_fun_Cm);
int_conv_Ia = sum(4/M*dalpha.*ind_fun_Ia);
int_conv_Iq = sum(2*b/(U*M)*dalphadot.*ind_fun_Iq);    % Notice that it should 1/M*dq correlates better than 2*b/(U*M)*dalphadot
int_conv_IM = sum(4*alpha/M^2*dM.*ind_fun_IM);
int_conv_Iam = sum(-1/M*dalpha.*ind_fun_Iam);
int_conv_Iqm = sum(-7*2*b/(12*a_inf*M^2)*dalphadot.*ind_fun_Iqm);
int_conv_IMm = sum(-alpha/M^2*dM.*ind_fun_IMm);

% Effective circulatory variables
w_tqc_eff = (cnaw_tqc-int_conv_C)/c_n_alpha;
Ucm_qs_eff = Ucm_qs-int_conv_Cm;

end