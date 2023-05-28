clc
clear
close all

%%
filename = 'gBest_016.mat';
load(filename,'best_struct');

disp("alpha_0L = " + num2str(rad2deg(best_struct.alpha_0L)));
disp("alpha_ds0 = " + num2str(rad2deg(best_struct.alpha_ds0)));
disp("alpha_ss = " + num2str(rad2deg(best_struct.alpha_ss)));
disp("alpha1_0 = " + num2str(rad2deg(best_struct.alpha1_0)));
disp("gamma_LS = " + num2str(best_struct.gamma_LS));
disp("delta_alpha_0 = " + num2str(rad2deg(best_struct.delta_alpha_0)));
disp("delta_alpha_1 = " + num2str(rad2deg(best_struct.delta_alpha_1)));
disp("delta_alpha_2 = " + num2str(best_struct.delta_alpha_2));
disp("nu_1 = " + num2str(best_struct.nu_1));
disp("nu_2 = " + num2str(best_struct.nu_2));
disp("c_d0 = " + num2str(best_struct.c_d0));
disp("c_m0 = " + num2str(best_struct.c_m0));
disp("c_n_alpha = " + num2str(deg2rad(best_struct.c_n_alpha)));
disp("d_cc = " + num2str(rad2deg(best_struct.d_cc)));
disp("d_cm = " + num2str(rad2deg(best_struct.d_cm)));
disp("E0 = " + num2str(best_struct.E0));
disp("E1 = " + num2str(best_struct.E1));
disp("g_v = " + num2str(best_struct.g_v));
disp("K0 = " + num2str(best_struct.K0));
disp("K1 = " + num2str(best_struct.K1));
disp("K2 = " + num2str(best_struct.K2));
disp("K3 = " + num2str(best_struct.K3));
disp("r0 = " + num2str(best_struct.r0));
disp("S1 = " + num2str(rad2deg(best_struct.S1)));
disp("S2 = " + num2str(rad2deg(best_struct.S2)));
disp("Ta = " + num2str(best_struct.Ta));
disp("Tf0 = " + num2str(best_struct.Tf0));
disp("TvL = " + num2str(best_struct.TvL));
disp("Vm = " + num2str(best_struct.Vm));
disp("Vn1 = " + num2str(best_struct.Vn1));
disp("Vn2 = " + num2str(best_struct.Vn2));
disp("z_cc = " + num2str(best_struct.z_cc));
disp("z_cm = " + num2str(best_struct.z_cm));
disp("gamma_TvL = " + num2str(best_struct.gamma_TvL));
disp("kappa = " + num2str(best_struct.kappa));
disp("df0_c = " + num2str(best_struct.df0_c));
disp("f0 = " + num2str(best_struct.f0));
disp("fb = " + num2str(best_struct.fb));
disp("fSig1n = " + num2str(best_struct.fSig1n));
disp("fSig1c = " + num2str(best_struct.fSig1c));
disp("fSig2n = " + num2str(best_struct.fSig2n));
disp("fS2n_ulpr = " + num2str(best_struct.fS2n_ulpr));
disp("fS2n_dlpr = " + num2str(best_struct.fS2n_dlpr));
disp("fS1n_u = " + num2str(best_struct.fS1n_u));
disp("fS1m_u = " + num2str(best_struct.fS1m_u));
disp("fS1c_u = " + num2str(best_struct.fS1c_u));
disp("fS1n_d = " + num2str(best_struct.fS1n_d));
disp("fS1m_d = " + num2str(best_struct.fS1m_d));
disp("fS1c_d = " + num2str(best_struct.fS1c_d));
disp("fS1n_ud = " + num2str(best_struct.fS1n_ud));
disp("fS1m_ud = " + num2str(best_struct.fS1m_ud));
disp("fS1c_ud = " + num2str(best_struct.fS1c_ud));
disp("fS2n_u = " + num2str(best_struct.fS2n_u));
disp("fS2m_u = " + num2str(best_struct.fS2m_u));
disp("fS2c_u = " + num2str(best_struct.fS2c_u));
disp("fS2n_d = " + num2str(best_struct.fS2n_d));
disp("fS2m_d = " + num2str(best_struct.fS2m_d));
disp("fS2c_d = " + num2str(best_struct.fS2c_d));
disp("fSS1 = " + num2str(best_struct.fSS1));
disp("fSS2 = " + num2str(best_struct.fSS2));
disp("K1_f = " + num2str(best_struct.K1_f));
