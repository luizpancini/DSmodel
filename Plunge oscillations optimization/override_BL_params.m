function [params,ps] = override_BL_params(p,params,model)

switch model
        case "BL"
            % Mach-dependent parameters
            ps.alpha_0L = p(1);
            ps.alpha_ds0 = p(2);
            ps.alpha_ss = p(3);
            ps.alpha1_0 = p(4);
            ps.gamma_LS = p(5);
            ps.delta_alpha_0 = p(6);
            ps.delta_alpha_1 = p(7);
            ps.delta_alpha_2 = p(8);
            ps.nu_1 = p(9);
            ps.nu_2 = p(10);
            ps.c_d0 = p(11);
            ps.c_m0 = p(12);
            ps.c_n_alpha = p(13);
            ps.d_cc = p(14);
            ps.d_cm = p(15);
            ps.E0 = p(16);
            ps.E1 = p(17);
            ps.g_v = p(18);
            ps.K0 = p(19);
            ps.K1 = p(20);
            ps.K2 = p(21);
            ps.K3 = p(22);
            ps.r0 = p(23);
            ps.S1 = p(24);
            ps.S2 = p(25);
            ps.Ta = p(26);
            ps.Tf0 = p(27);
            ps.TvL = p(28);
            ps.Vm = p(29);
            ps.Vn1 = p(30);
            ps.Vn2 = p(31);
            ps.z_cc = p(32);
            ps.z_cm = p(33);
            % Airfoil dependent
            ps.gamma_TvL = p(34);
            ps.kappa = p(35);
            ps.df0_c = p(36);
            ps.f0 = p(37);
            ps.fb = p(38);
            ps.fSig1n = p(39);
            ps.fSig1c = p(40);
            ps.fSig2n = p(41);
            ps.fS2n_ulpr = p(42);
            ps.fS2n_dlpr = p(43);
            ps.fS1n_u = p(44);
            ps.fS1m_u = p(45);
            ps.fS1c_u = p(46);
            ps.fS1n_d = p(47);
            ps.fS1m_d = p(48);
            ps.fS1c_d = p(49);
            ps.fS1n_ud = p(50);
            ps.fS1m_ud = p(51);
            ps.fS1c_ud = p(52);
            ps.fS2n_u = p(53);
            ps.fS2m_u = p(54);
            ps.fS2c_u = p(55);
            ps.fS2n_d = p(56);
            ps.fS2m_d = p(57);
            ps.fS2c_d = p(58);
            ps.fSS1 = p(59);
            ps.fSS2 = p(60);
            ps.K1_f = p(61);
            % Time delay constants adjustment for dimensional time
            ps.Ta = ps.Ta*params.b/params.U;
            ps.Tf0 = ps.Tf0*params.b/params.U;
            ps.TvL = ps.TvL*params.b/params.U;
            % Aerodynamic center
            ps.x_ac = 0.25-ps.K0;
            % Copy p struct to the params struct
            for fn = fieldnames(ps)'
                params.(fn{1}) = ps.(fn{1});
            end
            % Reset time delay constants to non-dimensional time for output
            ps.Ta = ps.Ta*params.U/params.b;
            ps.Tf0 = ps.Tf0*params.U/params.b;
            ps.TvL = ps.TvL*params.U/params.b;
end

end