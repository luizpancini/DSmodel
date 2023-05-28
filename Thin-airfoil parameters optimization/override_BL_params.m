function [params,ps] = override_BL_params(p,params,model)

switch model
        case "BL"
            % Mach-dependent parameters
            ps.alpha_0L =       p(1);
            ps.alpha_ds0 =      p(2);
            ps.alpha_ss =       p(3);
            ps.alpha1_0n =      p(4);
            ps.alpha1_0m =      p(5);
            ps.alpha1_0c =      p(6);
            ps.beta_Sig1n =     p(7);
            ps.beta_Sig1c =     p(8);
            ps.beta_Sig2n =     p(9);
            ps.beta_S2n_lpr =   p(10);
            ps.beta_S2c_lpr =   p(11);
            ps.beta_S1n_u =     p(12);
            ps.beta_S1m_u =     p(13);
            ps.beta_S1c_u =     p(14);
            ps.beta_S1n_d =     p(15);
            ps.beta_S1m_d =     p(16);
            ps.beta_S1c_d =     p(17);
            ps.beta_S2n_u =     p(18);
            ps.beta_S2m_u =     p(19);
            ps.beta_S2c_u =     p(20);
            ps.beta_S2n_d =     p(21);
            ps.beta_S2m_d =     p(22);
            ps.beta_S2c_d =     p(23);
            ps.gamma_LS =       p(24);
            ps.delta_alpha_0 =  p(25);
            ps.delta_alpha_1 =  p(26);
            ps.eta =            p(27);
            ps.kappa_0 =        p(28);
            ps.kappa_1 =        p(29);
            ps.kappa_2 =        p(30);
            ps.kappa_3 =        p(31);
            ps.lambda_1 =       p(32);
            ps.lambda_2 =       p(33);
            ps.mu_v2 =          p(34);
            ps.nu_1 =           p(35);
            ps.nu_2 =           p(36);
            ps.nu_3 =           p(37);
            ps.nu_4 =           p(38);
            ps.nu_5 =           p(39);
            ps.chi_u =          p(40);
            ps.chi_d =          p(41);
            ps.xi =             p(42);
            ps.zeta_a =         p(43);
            ps.c_d0 =           p(44);
            ps.c_m0 =           p(45);
            ps.c_n_alpha =      p(46);
            ps.d_cc =           p(47);
            ps.d_cm =           p(48);
            ps.E0 =             p(49);
            ps.E1 =             p(50);
            ps.f0_n =           p(51);
            ps.f0_m =           p(52);
            ps.f0_c =           p(53);
            ps.fb_n =           p(54);
            ps.fb_m =           p(55);
            ps.fb_c =           p(56);
            ps.g_v =            p(57);
            ps.g_v2 =           p(58);
            ps.K0 =             p(59);
            ps.K1 =             p(60);
            ps.K2 =             p(61);
            ps.r0 =             p(62);
            ps.S1_n =           p(63);
            ps.S1_m =           p(64);
            ps.S1_c =           p(65);
            ps.S2_n =           p(66);
            ps.S2_m =           p(67);
            ps.S2_c =           p(68);
            ps.Ta =             p(69);
            ps.Tf =             p(70);
            ps.Tv =             p(71);
            ps.Tv2 =            p(72);
            ps.Vm =             p(73);
            ps.Vc =             p(74);
            ps.Vn1 =            p(75);
            ps.Vn2 =            p(76);
            ps.Vn3 =            p(77);
            ps.z_cc =           p(78);
            ps.z_cm =           p(79);
            ps.A1 =             p(80);
            ps.A2 =           1-ps.A1;
            ps.b1 =             p(81);
            ps.b2 =             p(82);
            % Time delay constants adjustment for dimensional time
            ps.Ta = ps.Ta*params.b/params.U;
            ps.Tf = ps.Tf*params.b/params.U;
            ps.Tv = ps.Tv*params.b/params.U;
            ps.Tv2 = ps.Tv2*params.b/params.U;
            % Copy p struct to the params struct
            for fn = fieldnames(ps)'
                params.(fn{1}) = ps.(fn{1});
            end
            % Reset time delay constants to non-dimensional time for output
            ps.Ta = ps.Ta*params.U/params.b;
            ps.Tf = ps.Tf*params.U/params.b;
            ps.Tv = ps.Tv*params.U/params.b;
            ps.Tv2 = ps.Tv2*params.U/params.b;
end

end