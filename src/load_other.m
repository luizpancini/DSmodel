function [params,data] = load_other(case_now)

%% Load experimental data
switch case_now
    case 1
        % Test conditions
        M = 0.3;
        a_0 = 9.9*pi/180;
        a_1 = 9.9*pi/180;
        k = 0.1;
        a_1h = 0; k_h = 0;
        b = 0.61/2;             % Assumed, since they are comparing with data from McAlister (1982)
        a_inf = 340;            % Assumed MSL conditions
        ah = -1/2;              % Semichord-normalized pitch axis position after midchord
        airfoil = 'NACA0012';
        % McAlister (1982) [taken directly from Leishman and Crouse (1989)]
        file1 = '../Other Data/CN_alpha_LeishmanCrouse1989_A0_10_A1_10_EXP.mat';
        file2 = '../Other Data/CM_alpha_LeishmanCrouse1989_A0_10_A1_10_EXP.mat';
        load(file1)
        load(file2)
        alpha_exp_cn = CN_alpha_LeishmanCrouse1989_A0_10_A1_10_EXP(:,1);
        cn_exp = CN_alpha_LeishmanCrouse1989_A0_10_A1_10_EXP(:,2);
        alpha_exp_cm = CM_alpha_LeishmanCrouse1989_A0_10_A1_10_EXP(:,1);
        cm_exp = CM_alpha_LeishmanCrouse1989_A0_10_A1_10_EXP(:,2);
        alpha_exp_cl = nan; cl_exp = nan; alpha_exp_cd = nan; cd_exp = nan; alpha_exp_cc = nan; cc_exp = nan;
        authors{1} = 'McAlister et al. (1982)';
        % Leishman and Crouse (1989)
        file3 = '../Other Data/CN_alpha_LeishmanCrouse1989_A0_10_A1_10_BL.mat';
        file4 = '../Other Data/CM_alpha_LeishmanCrouse1989_A0_10_A1_10_BL.mat';
        load(file3)
        load(file4)
        alpha_mod_cn = CN_alpha_LeishmanCrouse1989_A0_10_A1_10_BL(:,1);
        cn_mod = CN_alpha_LeishmanCrouse1989_A0_10_A1_10_BL(:,2);
        alpha_mod_cm = CM_alpha_LeishmanCrouse1989_A0_10_A1_10_BL(:,1);
        cm_mod = CM_alpha_LeishmanCrouse1989_A0_10_A1_10_BL(:,2);
        alpha_mod_cl = nan; cl_mod = nan; alpha_mod_cd = nan; cd_mod = nan; alpha_mod_cc = nan; cc_mod = nan;
        authors{2} = 'Leishman and Crouse (1989)';
    case 2 % Actually A0 = 15 and A1 = 10 !!!
        % Test conditions
        M = 0.3;
        a_0 = 15*pi/180;
        a_1 = 10*pi/180;
        k = 0.1;
        a_1h = 0; k_h = 0;
        b = 0.61/2;             % Assumed, since they are comparing with data from McAlister (1982)
        a_inf = 340;            % Assumed MSL conditions
        ah = -1/2;              % Semichord-normalized pitch axis position after midchord
        airfoil = 'NACA0012';
        % McAllister (1982) [taken directly from Leishman and Crouse (1989)]
        file1 = '../Other Data/CN_alpha_LeishmanCrouse1989_A0_10_A1_15_EXP.mat';
        file2 = '../Other Data/CM_alpha_LeishmanCrouse1989_A0_10_A1_15_EXP.mat';
        load(file1)
        load(file2)
        alpha_exp_cn = CN_alpha_LeishmanCrouse1989_A0_10_A1_15_EXP(:,1);
        cn_exp = CN_alpha_LeishmanCrouse1989_A0_10_A1_15_EXP(:,2);
        alpha_exp_cm = CM_alpha_LeishmanCrouse1989_A0_10_A1_15_EXP(:,1);
        cm_exp = CM_alpha_LeishmanCrouse1989_A0_10_A1_15_EXP(:,2);
        alpha_exp_cl = nan; cl_exp = nan; alpha_exp_cd = nan; cd_exp = nan; alpha_exp_cc = nan; cc_exp = nan;
        authors{1} = 'McAlister et al. (1982)';
        % Leishman and Crouse (1989)
        file3 = '../Other Data/CN_alpha_LeishmanCrouse1989_A0_10_A1_15_BL.mat';
        file4 = '../Other Data/CM_alpha_LeishmanCrouse1989_A0_10_A1_15_BL.mat';
        load(file3)
        load(file4)
        alpha_mod_cn = CN_alpha_LeishmanCrouse1989_A0_10_A1_15_BL(:,1);
        cn_mod = CN_alpha_LeishmanCrouse1989_A0_10_A1_15_BL(:,2);
        alpha_mod_cm = CM_alpha_LeishmanCrouse1989_A0_10_A1_15_BL(:,1);
        cm_mod = CM_alpha_LeishmanCrouse1989_A0_10_A1_15_BL(:,2);
        alpha_mod_cl = nan; cl_mod = nan; alpha_mod_cd = nan; cd_mod = nan; alpha_mod_cc = nan; cc_mod = nan;
        authors{2} = 'Leishman and Crouse (1989)';
    case 3 % Sheng et al (2008) - A Modified Dynamic Stall Model for Low Mach Numbers - Fig. 9
        % Test conditions
        M = 0.12;
        a_0 = 14.7*pi/180;
        a_1 = 9.9*pi/180;
        k = 0.124;
        a_1h = 0; k_h = 0;
        authors{1} = 'EXP - Sheng et al. (2008)';
        authors{2} = 'MOD - Sheng et al. (2008)';
        b = 0.55/2;             % As displayed in Table 1 of Sheng et al (2006) - A New Stall-Onset Criterion for Low Speed Dynamic-Stall
        a_inf = 340;            % Assumed MSL conditions
        ah = -1/2;              % Semichord-normalized pitch axis position after midchord
        airfoil = 'NACA0012';
        % Experimental data
        file1 = '../Other Data/CN_alpha_Sheng_A0_15_A1_10_k0124_M0012_EXP.mat';
        file2 = '../Other Data/CM_alpha_Sheng_A0_15_A1_10_k0124_M0012_EXP.mat';
        file3 = '../Other Data/CC_alpha_Sheng_A0_15_A1_10_k0124_M0012_EXP.mat';
        file4 = '../Other Data/CD_alpha_Sheng_A0_15_A1_10_k0124_M0012_EXP.mat';
        load(file1)
        load(file2)
        load(file3)
        load(file4)
        alpha_exp_cn = CN_alpha_Sheng_A0_15_A1_10_k0124_M0012_EXP(:,1);
        cn_exp = CN_alpha_Sheng_A0_15_A1_10_k0124_M0012_EXP(:,2);
        alpha_exp_cm = CM_alpha_Sheng_A0_15_A1_10_k0124_M0012_EXP(:,1);
        cm_exp = CM_alpha_Sheng_A0_15_A1_10_k0124_M0012_EXP(:,2);
        alpha_exp_cd = CD_alpha_Sheng_A0_15_A1_10_k0124_M0012_EXP(:,1);
        cd_exp = CD_alpha_Sheng_A0_15_A1_10_k0124_M0012_EXP(:,2);
        alpha_exp_cc = CC_alpha_Sheng_A0_15_A1_10_k0124_M0012_EXP(:,1);
        cc_exp = CC_alpha_Sheng_A0_15_A1_10_k0124_M0012_EXP(:,2);
        alpha_exp_cl = nan; cl_exp = nan;
        % Model
        file1 = '../Other Data/CN_alpha_Sheng_A0_15_A1_10_k0124_M0012_MOD.mat';
        file2 = '../Other Data/CM_alpha_Sheng_A0_15_A1_10_k0124_M0012_MOD.mat';
        file3 = '../Other Data/CC_alpha_Sheng_A0_15_A1_10_k0124_M0012_MOD.mat';
        file4 = '../Other Data/CD_alpha_Sheng_A0_15_A1_10_k0124_M0012_MOD.mat';
        load(file1)
        load(file2)
        load(file3)
        load(file4)
        alpha_mod_cn = CN_alpha_Sheng_A0_15_A1_10_k0124_M0012_MOD(:,1);
        cn_mod = CN_alpha_Sheng_A0_15_A1_10_k0124_M0012_MOD(:,2);
        alpha_mod_cm = CM_alpha_Sheng_A0_15_A1_10_k0124_M0012_MOD(:,1);
        cm_mod = CM_alpha_Sheng_A0_15_A1_10_k0124_M0012_MOD(:,2);
        alpha_mod_cd = CD_alpha_Sheng_A0_15_A1_10_k0124_M0012_MOD(:,1);
        cd_mod = CD_alpha_Sheng_A0_15_A1_10_k0124_M0012_MOD(:,2);
        alpha_mod_cc = CC_alpha_Sheng_A0_15_A1_10_k0124_M0012_MOD(:,1);
        cc_mod = CC_alpha_Sheng_A0_15_A1_10_k0124_M0012_MOD(:,2);
        alpha_mod_cl = nan; cl_mod = nan;
    case 4 % Sheng et al (2007) - Improved Dynamic-Stall-Onset Criterion at Low Mach Numbers - Fig. 5
        % Test conditions
        M = 0.12;
        a_0 = 9.7*pi/180;
        a_1 = 9.9*pi/180;
        k = 0.025;
        a_1h = 0; k_h = 0;
        authors{1} = 'Sheng et al. (2007) - EXP';
        authors{2} = 'Sheng et al. (2007) - MOD';
        b = 0.55/2;             % As displayed in Table 1 of Sheng et al (2006) - A New Stall-Onset Criterion for Low Speed Dynamic-Stall
        a_inf = 340;            % Assumed MSL conditions
        ah = -1/2;              % Semichord-normalized pitch axis position after midchord
        airfoil = 'NACA0012';
        % Experimental
        file1 = '../Other Data/CN_alpha_Sheng_A0_10_A1_10_k0025_M0012_EXP.mat';
        load(file1)
        alpha_exp_cn = CN_alpha_Sheng_A0_10_A1_10_k0025_M0012_EXP(:,1);
        cn_exp = CN_alpha_Sheng_A0_10_A1_10_k0025_M0012_EXP(:,2);
        alpha_exp_cl = nan; cl_exp = nan; alpha_exp_cm = nan; cm_exp = nan; alpha_exp_cd = nan; cd_exp = nan; alpha_exp_cc = nan; cc_exp = nan;
        % Model
        file2 = '../Other Data/CN_alpha_Sheng_A0_10_A1_10_k0025_M0012_MOD.mat';
        load(file2)
        alpha_mod_cn = CN_alpha_Sheng_A0_10_A1_10_k0025_M0012_MOD(:,1);
        cn_mod = CN_alpha_Sheng_A0_10_A1_10_k0025_M0012_MOD(:,2);
        alpha_mod_cl = nan; cl_mod = nan; alpha_mod_cm = nan; cm_mod = nan; alpha_mod_cd = nan; cd_mod = nan; alpha_mod_cc = nan; cc_mod = nan;
    case 5 % By McCroskey and Pucci (1982) - Viscous-Inviscid interaction on Oscillating Airfoils - fig. 9
        % Test conditions
        M = 0.3;
        a_0 = 9*pi/180;
        a_1 = 5*pi/180;
        k = 0.20;
        a_1h = 0; k_h = 0;
        b = 0.61/2;             % Assumed
        a_inf = 340;            % Assumed MSL conditions
        ah = -1/2;              % Semichord-normalized pitch axis position after midchord
        airfoil = 'NACA0012';
        file1 = '../Other Data/NACA12_M030_alpha9p5_k020_CL.mat';
        file2 = '../Other Data/NACA12_M030_alpha9p5_k020_CM.mat';
        load(file1)
        load(file2)
        alpha_exp_cl = NACA12_M030_alpha9p5_k020_CL(:,1);
        cl_exp = NACA12_M030_alpha9p5_k020_CL(:,2);
        alpha_exp_cm = NACA12_M030_alpha9p5_k020_CM(:,1);
        cm_exp = NACA12_M030_alpha9p5_k020_CM(:,2);
        alpha_exp_cd = nan; cd_exp = nan; alpha_exp_cn = nan; cn_exp = nan; alpha_exp_cc = nan; cc_exp = nan;
        alpha_mod_cl = nan; cl_mod = nan; alpha_mod_cm = nan; cm_mod = nan; alpha_mod_cd = nan; cd_mod = nan; alpha_mod_cn = nan; cn_mod = nan; alpha_mod_cc = nan; cc_mod = nan;
        authors{1} = 'McCroskey and Pucci (1982)';
    case 6 % By McAlister et al - Dynamic Stall Experiments on the NACA 0012 Airfoil - (1978)
        % Test conditions
        M = 0.09;
        a_0 = 15*pi/180;
        a_1 = 10*pi/180;
        k = 0.15;
        a_1h = 0; k_h = 0;
        b = 1.22/2;
        a_inf = 340;            % Assumed MSL conditions
        ah = -1/2;              % Semichord-normalized pitch axis position after midchord
        airfoil = 'NACA0012';
        file1 = '../Other Data/CN_alpha_M_0090_k_0150_A0_15_A1_10.mat';
        file2 = '../Other Data/CM_alpha_M_0090_k_0150_A0_15_A1_10.mat';
        file3 = '../Other Data/CC_alpha_M_0090_k_0150_A0_15_A1_10.mat';
        load(file1)
        load(file2)
        load(file3)
        alpha_exp_cn = CN_alpha_M_0090_k_0150_A0_15_A1_10(:,1);
        cn_exp = CN_alpha_M_0090_k_0150_A0_15_A1_10(:,2);
        alpha_exp_cm = CM_alpha_M_0090_k_0150_A0_15_A1_10(:,1);
        cm_exp = CM_alpha_M_0090_k_0150_A0_15_A1_10(:,2);
        alpha_exp_cc = CC_alpha_M_0090_k_0150_A0_15_A1_10(:,1);
        cc_exp = -CC_alpha_M_0090_k_0150_A0_15_A1_10(:,2);
        alpha_exp_cl = nan; cl_exp = nan; alpha_exp_cd = nan; cd_exp = nan;
        alpha_mod_cl = nan; cl_mod = nan; alpha_mod_cm = nan; cm_mod = nan; alpha_mod_cd = nan; cd_mod = nan; alpha_mod_cn = nan; cn_mod = nan; alpha_mod_cc = nan; cc_mod = nan;
        authors{1} = 'McAlister et al. (1978)';
    case 7 % PARKER - Force and Pressure Measurements on an Airfoil Oscillating Trhough Stall - (1976)
        % Test conditions
        M = 0.075;
        a_0 = 15.6*pi/180;
        a_1 = 10.2*pi/180;
        k = 0.15;
        a_1h = 0; k_h = 0;
        authors{1} = 'Parker (1976) - closed';
        authors{2} = 'Parker (1976) - 2% open';
        b = 1.22/2;             % 4-foot chord
        a_inf = 340;            % Assumed MSL conditions
        ah = -1/2;              % Semichord-normalized pitch axis position after midchord
        airfoil = 'NACA0012';
        % Experimental data
        file1 = '../Other Data/CN_alpha_Re_2e6_k_0150_A0_16_A1_10_open.mat';
        file2 = '../Other Data/CN_alpha_Re_2e6_k_0150_A0_16_A1_10_closed.mat';
        file3 = '../Other Data/CM_alpha_Re_2e6_k_0150_A0_16_A1_10_open.mat';
        file4 = '../Other Data/CM_alpha_Re_2e6_k_0150_A0_16_A1_10_closed.mat';
        load(file1)
        load(file2)
        load(file3)
        load(file4)
        alpha_exp_cn = CN_alpha_Re_2e6_k_0150_A0_16_A1_10_closed(:,1);
        cn_exp = CN_alpha_Re_2e6_k_0150_A0_16_A1_10_closed(:,2);
        alpha_exp_cm = CM_alpha_Re_2e6_k_0150_A0_16_A1_10_closed(:,1);
        cm_exp = CM_alpha_Re_2e6_k_0150_A0_16_A1_10_closed(:,2);
        alpha_exp_cd = nan; cd_exp = nan; alpha_exp_cc = nan; cc_exp = nan; alpha_exp_cl = nan; cl_exp = nan;
        alpha_mod_cn = CN_alpha_Re_2e6_k_0150_A0_16_A1_10_open(:,1);
        cn_mod = CN_alpha_Re_2e6_k_0150_A0_16_A1_10_open(:,2);
        alpha_mod_cm = CM_alpha_Re_2e6_k_0150_A0_16_A1_10_open(:,1);
        cm_mod = CM_alpha_Re_2e6_k_0150_A0_16_A1_10_open(:,2);
        alpha_mod_cd = nan; cd_mod = nan; alpha_mod_cc = nan; cc_mod = nan;  alpha_mod_cl = nan; cl_mod = nan;
    case 8 % Sheng et al (2008) - A Modified Dynamic Stall Model for Low Mach Numbers - Fig. 10
        % Test conditions
        M = 0.12;
        a_0 = 15.3*pi/180;
        a_1 = 9.9*pi/180;
        k = 0.074;
        a_1h = 0; k_h = 0;
        authors{1} = 'EXP - Sheng et al. (2008)';
        authors{2} = 'MOD - Sheng et al. (2008)';
        b = 0.55/2;             % As displayed in Table 1 of Sheng et al (2006) - A New Stall-Onset Criterion for Low Speed Dynamic-Stall
        a_inf = 340;            % Assumed MSL conditions
        ah = -1/2;              % Semichord-normalized pitch axis position after midchord
        airfoil = 'S809';
        % Experimental data
        file1 = '../Other Data/CN_alpha_Sheng_S809_A0_15_A1_10_k0074_M0012_EXP.mat';
        file2 = '../Other Data/CM_alpha_Sheng_S809_A0_15_A1_10_k0074_M0012_EXP.mat';
        file3 = '../Other Data/CC_alpha_Sheng_S809_A0_15_A1_10_k0074_M0012_EXP.mat';
        file4 = '../Other Data/CD_alpha_Sheng_S809_A0_15_A1_10_k0074_M0012_EXP.mat';
        load(file1)
        load(file2)
        load(file3)
        load(file4)
        alpha_exp_cn = CN_alpha_Sheng_S809_A0_15_A1_10_k0074_M0012_EXP(:,1);
        cn_exp = CN_alpha_Sheng_S809_A0_15_A1_10_k0074_M0012_EXP(:,2);
        alpha_exp_cm = CM_alpha_Sheng_S809_A0_15_A1_10_k0074_M0012_EXP(:,1);
        cm_exp = CM_alpha_Sheng_S809_A0_15_A1_10_k0074_M0012_EXP(:,2);
        alpha_exp_cd = CD_alpha_Sheng_S809_A0_15_A1_10_k0074_M0012_EXP(:,1);
        cd_exp = CD_alpha_Sheng_S809_A0_15_A1_10_k0074_M0012_EXP(:,2);
        alpha_exp_cc = CC_alpha_Sheng_S809_A0_15_A1_10_k0074_M0012_EXP(:,1);
        cc_exp = CC_alpha_Sheng_S809_A0_15_A1_10_k0074_M0012_EXP(:,2);
        alpha_exp_cl = nan; cl_exp = nan;
        % Model
        file1 = '../Other Data/CN_alpha_Sheng_S809_A0_15_A1_10_k0074_M0012_MOD.mat';
        file2 = '../Other Data/CM_alpha_Sheng_S809_A0_15_A1_10_k0074_M0012_MOD.mat';
        file3 = '../Other Data/CC_alpha_Sheng_S809_A0_15_A1_10_k0074_M0012_MOD.mat';
        file4 = '../Other Data/CD_alpha_Sheng_S809_A0_15_A1_10_k0074_M0012_MOD.mat';
        load(file1)
        load(file2)
        load(file3)
        load(file4)
        alpha_mod_cn = CN_alpha_Sheng_S809_A0_15_A1_10_k0074_M0012_MOD(:,1);
        cn_mod = CN_alpha_Sheng_S809_A0_15_A1_10_k0074_M0012_MOD(:,2);
        alpha_mod_cm = CM_alpha_Sheng_S809_A0_15_A1_10_k0074_M0012_MOD(:,1);
        cm_mod = CM_alpha_Sheng_S809_A0_15_A1_10_k0074_M0012_MOD(:,2);
        alpha_mod_cd = CD_alpha_Sheng_S809_A0_15_A1_10_k0074_M0012_MOD(:,1);
        cd_mod = CD_alpha_Sheng_S809_A0_15_A1_10_k0074_M0012_MOD(:,2);
        alpha_mod_cc = CC_alpha_Sheng_S809_A0_15_A1_10_k0074_M0012_MOD(:,1);
        cc_mod = CC_alpha_Sheng_S809_A0_15_A1_10_k0074_M0012_MOD(:,2);
        alpha_mod_cl = nan; cl_mod = nan;
    case 9
        % McAlister (1982)
        [params,data] = load_NASA(14210);
        [M,U,b,a_inf,beta,a_0,a_1,k,airfoil,ah] = struct2vars(params);
        [authors,alpha_exp_cl,alpha_exp_cm,alpha_exp_cd,cl_exp,cm_exp,cd_exp] = struct2vars(data);
        a_1h = 0; k_h = 0;
        % SHUM & LEE (2020)
        file1 = '../Other Data/SHUM_LEE_CL.mat';
        file2 = '../Other Data/SHUM_LEE_CM.mat';
        file3 = '../Other Data/SHUM_LEE_CD.mat';
        load(file1)
        load(file2)
        load(file3)
        alpha_mod_cl = SHUM_LEE_CL(:,1);
        cl_mod = SHUM_LEE_CL(:,2);
        alpha_mod_cm = SHUM_LEE_CM(:,1);
        cm_mod = SHUM_LEE_CM(:,2);
        alpha_mod_cd = SHUM_LEE_CD(:,1);
        cd_mod = SHUM_LEE_CD(:,2);
        alpha_exp_cn = nan; cn_exp = nan; alpha_exp_cc = nan; cc_exp = nan;
        alpha_mod_cn = nan; cn_mod = nan; alpha_mod_cc = nan; cc_mod = nan;
        authors{2} = 'Shum & Lee (2020)';
    case 10
        % Test conditions
        M = 0.075; 
        a_0 = 17.5*pi/180;
        a_1 = 22.5*pi/180;
        k = 1e-3;
        a_1h = 0; k_h = 0;
        authors{1} = 'Ramsay et al. (1996)';
        authors{1} = 'NREL';
        b = 0.437/2;
        a_inf = 340;
        ah = -1/2;
        airfoil = 'S809';
        % Experimental data
        file = '../Other Data/S809_static_M_0075.mat';
        load(file,'alpha','c_l','c_d','c_m','c_n','c_c')
        alpha_exp_cl = alpha; alpha_exp_cd = alpha; alpha_exp_cn = alpha; alpha_exp_cm = alpha; alpha_exp_cc = alpha;
        cl_exp = c_l; cd_exp = c_d; cn_exp = c_n; cm_exp = c_m; cc_exp = c_c;
        alpha_mod_cl = nan; cl_mod = nan; alpha_mod_cm = nan; cm_mod = nan; alpha_mod_cd = nan; cd_mod = nan; alpha_mod_cn = nan; cn_mod = nan; alpha_mod_cc = nan; cc_mod = nan;
    case 11
        % Test conditions
        M = 0.150; 
        a_0 = 10.0*pi/180;
        a_1 = 15.0*pi/180;
        k = 5e-4;
        a_1h = 0; k_h = 0;
        authors{1} = 'Ramsay et al. (1996)';
        authors{1} = 'NREL';
        b = 0.437/2;
        a_inf = 340;
        ah = -1/2;
        airfoil = 'S809';
        % Experimental data
        file = '../Other Data/S809_static_M_0150.mat';
        load(file,'alpha','c_l','c_d','c_m','c_n','c_c')
        alpha_exp_cl = alpha; alpha_exp_cd = alpha; alpha_exp_cn = alpha; alpha_exp_cm = alpha; alpha_exp_cc = alpha;
        cl_exp = c_l; cd_exp = c_d; cn_exp = c_n; cm_exp = c_m; cc_exp = c_c;
        alpha_mod_cl = nan; cl_mod = nan; alpha_mod_cm = nan; cm_mod = nan; alpha_mod_cd = nan; cd_mod = nan; alpha_mod_cn = nan; cn_mod = nan; alpha_mod_cc = nan; cc_mod = nan;    
    case 12
        % Test conditions
        M = 0.2; 
        a_0 = 12.5*pi/180;
        a_1 = 12.5*pi/180;
        k = 0.001;
        authors{1} = 'Liiva et al. (1968)';
        b = 6.38*0.0254/2;
        a_inf = 325;
        ah = -1/2;
        airfoil = 'Vertol-23010';
        a_1h = 0; k_h = 0;
        % Experimental data
        file1 = '../Other Data/CN_Vertol23010_steady_M_02.mat';
        file2 = '../Other Data/CM_Vertol23010_steady_M_02.mat';
        load(file1)
        load(file2)
        alpha_exp_cn = CN_Vertol23010_steady_M_02(:,1);
        cn_exp = CN_Vertol23010_steady_M_02(:,2);
        alpha_exp_cm = CM_Vertol23010_steady_M_02(:,1);
        cm_exp = CM_Vertol23010_steady_M_02(:,2);
        alpha_exp_cl = nan; cl_exp = nan; alpha_exp_cd = nan; cd_exp = nan; alpha_exp_cc = nan; cc_exp = nan;
        alpha_mod_cl = nan; cl_mod = nan; alpha_mod_cm = nan; cm_mod = nan; alpha_mod_cd = nan; cd_mod = nan; alpha_mod_cn = nan; cn_mod = nan; alpha_mod_cc = nan; cc_mod = nan;
    case 13
        % Test conditions
        M = 0.4;
        a_0 = 12.5*pi/180;
        a_1 = 12.5*pi/180;
        k = 0.001;
        authors{1} = 'Liiva et al. (1968)';
        b = 6.38*0.0254/2;
        a_inf = 325;
        ah = -1/2;
        airfoil = 'Vertol-23010';
        a_1h = 0; k_h = 0;
        % Experimental data
        file1 = '../Other Data/CN_Vertol23010_steady_M_04.mat';
        file2 = '../Other Data/CM_Vertol23010_steady_M_04.mat';
        load(file1)
        load(file2)
        alpha_exp_cn = CN_Vertol23010_steady_M_04(:,1);
        cn_exp = CN_Vertol23010_steady_M_04(:,2);
        alpha_exp_cm = CM_Vertol23010_steady_M_04(:,1);
        cm_exp = CM_Vertol23010_steady_M_04(:,2);
        alpha_exp_cl = nan; cl_exp = nan; alpha_exp_cd = nan; cd_exp = nan; alpha_exp_cc = nan; cc_exp = nan;
        alpha_mod_cl = nan; cl_mod = nan; alpha_mod_cm = nan; cm_mod = nan; alpha_mod_cd = nan; cd_mod = nan; alpha_mod_cn = nan; cn_mod = nan; alpha_mod_cc = nan; cc_mod = nan;
    case 14
        % Test conditions
        M = 0.4;
        a_0 = 0.05*pi/180;
        a_1 = 5.05*pi/180;
        k = 0.125;
        authors = {'Liiva et al. (1968)'; 'Tyler & Leishman (1992)'};
        b = 6.38*0.0254/2;
        a_inf = 325;
        ah = -1/2;
        airfoil = 'Vertol-23010';
        a_1h = 0; k_h = 0;
        % Experimental data
        file1 = '../Other Data/CN_Vertol23010_M_04_k_0125_a0_0_a1_5_exp.mat';
        file2 = '../Other Data/CN_Vertol23010_M_04_k_0125_a0_0_a1_5_mod.mat';
        file3 = '../Other Data/CM_Vertol23010_M_04_k_0125_a0_0_a1_5_exp.mat';
        file4 = '../Other Data/CM_Vertol23010_M_04_k_0125_a0_0_a1_5_mod.mat';
        load(file1)
        load(file2)
        load(file3)
        load(file4)
        alpha_exp_cn = CN_Vertol23010_M_04_k_0125_a0_0_a1_5_exp(:,1);
        cn_exp = CN_Vertol23010_M_04_k_0125_a0_0_a1_5_exp(:,2);
        alpha_exp_cm = CM_Vertol23010_M_04_k_0125_a0_0_a1_5_exp(:,1);
        cm_exp = CM_Vertol23010_M_04_k_0125_a0_0_a1_5_exp(:,2);
        alpha_mod_cn = CN_Vertol23010_M_04_k_0125_a0_0_a1_5_mod(:,1);
        cn_mod = CN_Vertol23010_M_04_k_0125_a0_0_a1_5_mod(:,2);
        alpha_mod_cm = CM_Vertol23010_M_04_k_0125_a0_0_a1_5_mod(:,1);
        cm_mod = CM_Vertol23010_M_04_k_0125_a0_0_a1_5_mod(:,2);
        alpha_exp_cl = nan; cl_exp = nan; alpha_exp_cd = nan; cd_exp = nan; alpha_exp_cc = nan; cc_exp = nan;
        alpha_mod_cl = nan; cl_mod = nan; alpha_mod_cd = nan; cd_mod = nan; alpha_mod_cc = nan; cc_mod = nan;
    case 15
        % Test conditions
        M = 0.4;
        a_0 = 12.22*pi/180;
        a_1 = 4.83*pi/180;
        k = 0.062;
        authors = {'Liiva et al. (1968)'; 'Tyler & Leishman (1992)'};
        b = 6.38*0.0254/2;
        a_inf = 325;
        ah = -1/2;
        airfoil = 'Vertol-23010';
        a_1h = 0; k_h = 0;
        % Experimental data
        file1 = '../Other Data/CN_Vertol23010_M_04_k_0062_a0_12_a1_5_exp.mat';
        file2 = '../Other Data/CN_Vertol23010_M_04_k_0062_a0_12_a1_5_mod.mat';
        file3 = '../Other Data/CM_Vertol23010_M_04_k_0062_a0_12_a1_5_exp.mat';
        file4 = '../Other Data/CM_Vertol23010_M_04_k_0062_a0_12_a1_5_mod.mat';
        load(file1)
        load(file2)
        load(file3)
        load(file4)
        alpha_exp_cn = CN_Vertol23010_M_04_k_0062_a0_12_a1_5_exp(:,1);
        cn_exp = CN_Vertol23010_M_04_k_0062_a0_12_a1_5_exp(:,2);
        alpha_exp_cm = CM_Vertol23010_M_04_k_0062_a0_12_a1_5_exp(:,1);
        cm_exp = CM_Vertol23010_M_04_k_0062_a0_12_a1_5_exp(:,2);
        alpha_mod_cn = CN_Vertol23010_M_04_k_0062_a0_12_a1_5_mod(:,1);
        cn_mod = CN_Vertol23010_M_04_k_0062_a0_12_a1_5_mod(:,2);
        alpha_mod_cm = CM_Vertol23010_M_04_k_0062_a0_12_a1_5_mod(:,1);
        cm_mod = CM_Vertol23010_M_04_k_0062_a0_12_a1_5_mod(:,2);
        alpha_exp_cl = nan; cl_exp = nan; alpha_exp_cd = nan; cd_exp = nan; alpha_exp_cc = nan; cc_exp = nan;
        alpha_mod_cl = nan; cl_mod = nan; alpha_mod_cd = nan; cd_mod = nan; alpha_mod_cc = nan; cc_mod = nan;     
   case 16
        % Test conditions
        M = 0.4;
        a_0 = 12.29*pi/180;
        a_1 = 4.93*pi/180;
        k = 0.124;
        authors = {'Liiva et al. (1968)'; 'Tyler & Leishman (1992)'};
        b = 6.38*0.0254/2;
        a_inf = 325;
        ah = -1/2;
        airfoil = 'Vertol-23010';
        a_1h = 0; k_h = 0;
        % Experimental data
        file1 = '../Other Data/CN_Vertol23010_M_04_k_0124_a0_12_a1_5_exp.mat';
        file2 = '../Other Data/CN_Vertol23010_M_04_k_0124_a0_12_a1_5_mod.mat';
        file3 = '../Other Data/CM_Vertol23010_M_04_k_0124_a0_12_a1_5_exp.mat';
        file4 = '../Other Data/CM_Vertol23010_M_04_k_0124_a0_12_a1_5_mod.mat';
        load(file1)
        load(file2)
        load(file3)
        load(file4)
        alpha_exp_cn = CN_Vertol23010_M_04_k_0124_a0_12_a1_5_exp(:,1);
        cn_exp = CN_Vertol23010_M_04_k_0124_a0_12_a1_5_exp(:,2);
        alpha_exp_cm = CM_Vertol23010_M_04_k_0124_a0_12_a1_5_exp(:,1);
        cm_exp = CM_Vertol23010_M_04_k_0124_a0_12_a1_5_exp(:,2);
        alpha_mod_cn = CN_Vertol23010_M_04_k_0124_a0_12_a1_5_mod(:,1);
        cn_mod = CN_Vertol23010_M_04_k_0124_a0_12_a1_5_mod(:,2);
        alpha_mod_cm = CM_Vertol23010_M_04_k_0124_a0_12_a1_5_mod(:,1);
        cm_mod = CM_Vertol23010_M_04_k_0124_a0_12_a1_5_mod(:,2);
        alpha_exp_cl = nan; cl_exp = nan; alpha_exp_cd = nan; cd_exp = nan; alpha_exp_cc = nan; cc_exp = nan;
        alpha_mod_cl = nan; cl_mod = nan; alpha_mod_cd = nan; cd_mod = nan; alpha_mod_cc = nan; cc_mod = nan;
    case 17
        % Test conditions
        M = 0.4;
        a_0 = 12.4*pi/180;
        a_1 = 5.2*pi/180;
        k = 0.185;
        authors = {'Liiva et al. (1968)'};
        b = 6.38*0.0254/2;
        a_inf = 325;
        ah = -1/2;
        airfoil = 'Vertol-23010';
        a_1h = 0; k_h = 0;
        % Experimental data
        file1 = '../Other Data/CN_Vertol23010_M_04_k_0185_a0_123_a1_5_exp.mat';
        file2 = '../Other Data/CM_Vertol23010_M_04_k_0185_a0_123_a1_5_exp.mat';
        load(file1)
        load(file2)
        alpha_exp_cn = CN_Vertol23010_M_04_k_0185_a0_123_a1_5_exp(:,1);
        cn_exp = CN_Vertol23010_M_04_k_0185_a0_123_a1_5_exp(:,2);
        alpha_exp_cm = CM_Vertol23010_M_04_k_0185_a0_123_a1_5_exp(:,1);
        cm_exp = CM_Vertol23010_M_04_k_0185_a0_123_a1_5_exp(:,2);
        alpha_mod_cn = nan; cn_mod = nan; alpha_mod_cm = nan; cm_mod = nan;
        alpha_exp_cl = nan; cl_exp = nan; alpha_exp_cd = nan; cd_exp = nan; alpha_exp_cc = nan; cc_exp = nan;
        alpha_mod_cl = nan; cl_mod = nan; alpha_mod_cd = nan; cd_mod = nan; alpha_mod_cc = nan; cc_mod = nan;    
    case 18
        % Test conditions
        M = 0.4;
        a_0 = 12.5*pi/180;
        a_1 = 5.3*pi/180;
        k = 0.238;
        authors = {'Liiva et al. (1968)'};
        b = 6.38*0.0254/2;
        a_inf = 325;
        ah = -1/2;
        airfoil = 'Vertol-23010';
        a_1h = 0; k_h = 0;
        % Experimental data
        file1 = '../Other Data/CN_Vertol23010_M_04_k_0238_a0_123_a1_5_exp.mat';
        file2 = '../Other Data/CM_Vertol23010_M_04_k_0238_a0_123_a1_5_exp.mat';
        load(file1)
        load(file2)
        alpha_exp_cn = CN_Vertol23010_M_04_k_0238_a0_123_a1_5_exp(:,1);
        cn_exp = CN_Vertol23010_M_04_k_0238_a0_123_a1_5_exp(:,2);
        alpha_exp_cm = CM_Vertol23010_M_04_k_0238_a0_123_a1_5_exp(:,1);
        cm_exp = CM_Vertol23010_M_04_k_0238_a0_123_a1_5_exp(:,2);
        alpha_mod_cn = nan; cn_mod = nan; alpha_mod_cm = nan; cm_mod = nan;
        alpha_exp_cl = nan; cl_exp = nan; alpha_exp_cd = nan; cd_exp = nan; alpha_exp_cc = nan; cc_exp = nan;
        alpha_mod_cl = nan; cl_mod = nan; alpha_mod_cd = nan; cd_mod = nan; alpha_mod_cc = nan; cc_mod = nan;  
    case 19
        % Test conditions
        M = 0.4;
        a_0 = 12.3*pi/180;
        a_1 = 5.6*pi/180;
        k = 0.302;
        authors = {'Liiva et al. (1968)'};
        b = 6.38*0.0254/2;
        a_inf = 325;
        ah = -1/2;
        airfoil = 'Vertol-23010';
        a_1h = 0; k_h = 0;
        % Experimental data
        file1 = '../Other Data/CN_Vertol23010_M_04_k_0302_a0_123_a1_5_exp.mat';
        file2 = '../Other Data/CM_Vertol23010_M_04_k_0302_a0_123_a1_5_exp.mat';
        load(file1)
        load(file2)
        alpha_exp_cn = CN_Vertol23010_M_04_k_0302_a0_123_a1_5_exp(:,1);
        cn_exp = CN_Vertol23010_M_04_k_0302_a0_123_a1_5_exp(:,2);
        alpha_exp_cm = CM_Vertol23010_M_04_k_0302_a0_123_a1_5_exp(:,1);
        cm_exp = CM_Vertol23010_M_04_k_0302_a0_123_a1_5_exp(:,2);
        alpha_mod_cn = nan; cn_mod = nan; alpha_mod_cm = nan; cm_mod = nan;
        alpha_exp_cl = nan; cl_exp = nan; alpha_exp_cd = nan; cd_exp = nan; alpha_exp_cc = nan; cc_exp = nan;
        alpha_mod_cl = nan; cl_mod = nan; alpha_mod_cd = nan; cd_mod = nan; alpha_mod_cc = nan; cc_mod = nan; 
    case 20
        % Test conditions
        M = 0.4;
        a_0 = 12.5*pi/180;
        a_1 = 5.65*pi/180;
        k = 0.355;
        authors = {'Liiva et al. (1968)'};
        b = 6.38*0.0254/2;
        a_inf = 325;
        ah = -1/2;
        airfoil = 'Vertol-23010';
        a_1h = 0; k_h = 0;
        % Experimental data
        file1 = '../Other Data/CN_Vertol23010_M_04_k_0355_a0_123_a1_5_exp.mat';
        file2 = '../Other Data/CM_Vertol23010_M_04_k_0355_a0_123_a1_5_exp.mat';
        load(file1)
        load(file2)
        alpha_exp_cn = CN_Vertol23010_M_04_k_0355_a0_123_a1_5_exp(:,1);
        cn_exp = CN_Vertol23010_M_04_k_0355_a0_123_a1_5_exp(:,2);
        alpha_exp_cm = CM_Vertol23010_M_04_k_0355_a0_123_a1_5_exp(:,1);
        cm_exp = CM_Vertol23010_M_04_k_0355_a0_123_a1_5_exp(:,2);
        alpha_mod_cn = nan; cn_mod = nan; alpha_mod_cm = nan; cm_mod = nan;
        alpha_exp_cl = nan; cl_exp = nan; alpha_exp_cd = nan; cd_exp = nan; alpha_exp_cc = nan; cc_exp = nan;
        alpha_mod_cl = nan; cl_mod = nan; alpha_mod_cd = nan; cd_mod = nan; alpha_mod_cc = nan; cc_mod = nan;
    case 21
        % Test conditions
        M = 0.4;
        a_0 = 15.07*pi/180;
        a_1 = 4.99*pi/180;
        k = 0.122;
        authors = {'Liiva et al. (1968)'; 'Tyler & Leishman (1992)'};
        b = 6.38*0.0254/2;
        a_inf = 325;
        ah = -1/2;
        airfoil = 'Vertol-23010';
        a_1h = 0; k_h = 0;
        % Experimental data
        file1 = '../Other Data/CN_Vertol23010_M_04_k_0122_a0_15_a1_5_exp.mat';
        file2 = '../Other Data/CN_Vertol23010_M_04_k_0122_a0_15_a1_5_mod.mat';
        file3 = '../Other Data/CM_Vertol23010_M_04_k_0122_a0_15_a1_5_exp.mat';
        file4 = '../Other Data/CM_Vertol23010_M_04_k_0122_a0_15_a1_5_mod.mat';
        load(file1)
        load(file2)
        load(file3)
        load(file4)
        alpha_exp_cn = CN_Vertol23010_M_04_k_0122_a0_15_a1_5_exp(:,1);
        cn_exp = CN_Vertol23010_M_04_k_0122_a0_15_a1_5_exp(:,2);
        alpha_exp_cm = CM_Vertol23010_M_04_k_0122_a0_15_a1_5_exp(:,1);
        cm_exp = CM_Vertol23010_M_04_k_0122_a0_15_a1_5_exp(:,2);
        alpha_mod_cn = CN_Vertol23010_M_04_k_0122_a0_15_a1_5_mod(:,1);
        cn_mod = CN_Vertol23010_M_04_k_0122_a0_15_a1_5_mod(:,2);
        alpha_mod_cm = CM_Vertol23010_M_04_k_0122_a0_15_a1_5_mod(:,1);
        cm_mod = CM_Vertol23010_M_04_k_0122_a0_15_a1_5_mod(:,2);
        alpha_exp_cl = nan; cl_exp = nan; alpha_exp_cd = nan; cd_exp = nan; alpha_exp_cc = nan; cc_exp = nan;
        alpha_mod_cl = nan; cl_mod = nan; alpha_mod_cd = nan; cd_mod = nan; alpha_mod_cc = nan; cc_mod = nan;
    case 22
        % Test conditions
        M = 0.4;
        a_0 = 9.85*pi/180;
        a_1 = 7.67*pi/180;
        k = 0.124;
        authors = {'Liiva et al. (1968)'};
        b = 6.38*0.0254/2;
        a_inf = 325;
        ah = -1/2;
        airfoil = 'Vertol-23010';
        a_1h = 0; k_h = 0;
        % Experimental data
        file1 = '../Other Data/CN_Vertol23010_M_04_k_0125_a0_10_a1_75_exp.mat';
        file2 = '../Other Data/CM_Vertol23010_M_04_k_0125_a0_10_a1_75_exp.mat';
        load(file1)
        load(file2)
        alpha_exp_cn = CN_Vertol23010_M_04_k_0125_a0_10_a1_75_exp(:,1);
        cn_exp = CN_Vertol23010_M_04_k_0125_a0_10_a1_75_exp(:,2);
        alpha_exp_cm = CM_Vertol23010_M_04_k_0125_a0_10_a1_75_exp(:,1);
        cm_exp = CM_Vertol23010_M_04_k_0125_a0_10_a1_75_exp(:,2);
        alpha_mod_cn = nan; cn_mod = nan; alpha_mod_cm = nan; cm_mod = nan;
        alpha_exp_cl = nan; cl_exp = nan; alpha_exp_cd = nan; cd_exp = nan; alpha_exp_cc = nan; cc_exp = nan;
        alpha_mod_cl = nan; cl_mod = nan; alpha_mod_cd = nan; cd_mod = nan; alpha_mod_cc = nan; cc_mod = nan;    
    case 23
        % Test conditions
        M = 0.4;
        a_0 = 14.6*pi/180;
        a_1 = 2.46*pi/180;
        k = 0.125;
        authors = {'Liiva et al. (1968)'};
        b = 6.38*0.0254/2;
        a_inf = 325;
        ah = -1/2;
        airfoil = 'Vertol-23010';
        a_1h = 0; k_h = 0;
        % Experimental data
        file1 = '../Other Data/CN_Vertol23010_M_04_k_0125_a0_15_a1_25_exp.mat';
        file2 = '../Other Data/CM_Vertol23010_M_04_k_0125_a0_15_a1_25_exp.mat';
        load(file1)
        load(file2)
        alpha_exp_cn = CN_Vertol23010_M_04_k_0125_a0_15_a1_25_exp(:,1);
        cn_exp = CN_Vertol23010_M_04_k_0125_a0_15_a1_25_exp(:,2);
        alpha_exp_cm = CM_Vertol23010_M_04_k_0125_a0_15_a1_25_exp(:,1);
        cm_exp = CM_Vertol23010_M_04_k_0125_a0_15_a1_25_exp(:,2);
        alpha_mod_cn = nan; cn_mod = nan; alpha_mod_cm = nan; cm_mod = nan;
        alpha_exp_cl = nan; cl_exp = nan; alpha_exp_cd = nan; cd_exp = nan; alpha_exp_cc = nan; cc_exp = nan;
        alpha_mod_cl = nan; cl_mod = nan; alpha_mod_cd = nan; cd_mod = nan; alpha_mod_cc = nan; cc_mod = nan;        
    case 24
        % Test conditions
        M = 0.4;
        a_0 = 0.26*pi/180;
        a_1h = 3.10*pi/180;
        k_h = 0.129;
        authors = {'Liiva et al. (1968)'; 'Tyler & Leishman (1992)'};
        b = 6.38*0.0254/2;
        a_inf = 325;
        ah = -1/2;
        airfoil = 'Vertol-23010';
        a_1 = 0; k = 0;
        % Experimental data
        file1 = '../Other Data/CN_Vertol23010_M_04_k_0129_a0_0_a1h_3_exp.mat';
        file2 = '../Other Data/CN_Vertol23010_M_04_k_0129_a0_0_a1h_3_mod.mat';
        file3 = '../Other Data/CM_Vertol23010_M_04_k_0129_a0_0_a1h_3_exp.mat';
        file4 = '../Other Data/CM_Vertol23010_M_04_k_0129_a0_0_a1h_3_mod.mat';
        load(file1)
        load(file2)
        load(file3)
        load(file4)
        alpha_exp_cn = CN_Vertol23010_M_04_k_0129_a0_0_a1h_3_exp(:,1);
        cn_exp = CN_Vertol23010_M_04_k_0129_a0_0_a1h_3_exp(:,2);
        alpha_exp_cm = CM_Vertol23010_M_04_k_0129_a0_0_a1h_3_exp(:,1);
        cm_exp = CM_Vertol23010_M_04_k_0129_a0_0_a1h_3_exp(:,2);
        alpha_mod_cn = CN_Vertol23010_M_04_k_0129_a0_0_a1h_3_mod(:,1);
        cn_mod = CN_Vertol23010_M_04_k_0129_a0_0_a1h_3_mod(:,2);
        alpha_mod_cm = CM_Vertol23010_M_04_k_0129_a0_0_a1h_3_mod(:,1);
        cm_mod = CM_Vertol23010_M_04_k_0129_a0_0_a1h_3_mod(:,2);
        alpha_exp_cl = nan; cl_exp = nan; alpha_exp_cd = nan; cd_exp = nan; alpha_exp_cc = nan; cc_exp = nan;
        alpha_mod_cl = nan; cl_mod = nan; alpha_mod_cd = nan; cd_mod = nan; alpha_mod_cc = nan; cc_mod = nan;
    case 25
        % Test conditions
        M = 0.4;
        a_0 = 12.53*pi/180;
        a_1h = 2.05*pi/180;
        k_h = 0.068; % Liiva's data shows 0.068 instead of Tyler's 0.058
        authors = {'Liiva et al. (1968)'; 'Tyler & Leishman (1992)'};
        b = 6.38*0.0254/2;
        a_inf = 325;
        ah = -1/2;
        airfoil = 'Vertol-23010';
        a_1 = 0; k = 0;
        % Experimental data
        file1 = '../Other Data/CN_Vertol23010_M_04_k_0058_a0_13_a1h_2_exp.mat';
        file2 = '../Other Data/CN_Vertol23010_M_04_k_0058_a0_13_a1h_2_mod.mat';
        file3 = '../Other Data/CM_Vertol23010_M_04_k_0058_a0_13_a1h_2_exp.mat';
        file4 = '../Other Data/CM_Vertol23010_M_04_k_0058_a0_13_a1h_2_mod.mat';
        load(file1)
        load(file2)
        load(file3)
        load(file4)
        alpha_exp_cn = CN_Vertol23010_M_04_k_0058_a0_13_a1h_2_exp(:,1);
        cn_exp = CN_Vertol23010_M_04_k_0058_a0_13_a1h_2_exp(:,2);
        alpha_exp_cm = CM_Vertol23010_M_04_k_0058_a0_13_a1h_2_exp(:,1);
        cm_exp = CM_Vertol23010_M_04_k_0058_a0_13_a1h_2_exp(:,2);
        alpha_mod_cn = CN_Vertol23010_M_04_k_0058_a0_13_a1h_2_mod(:,1);
        cn_mod = CN_Vertol23010_M_04_k_0058_a0_13_a1h_2_mod(:,2);
        alpha_mod_cm = CM_Vertol23010_M_04_k_0058_a0_13_a1h_2_mod(:,1);
        cm_mod = CM_Vertol23010_M_04_k_0058_a0_13_a1h_2_mod(:,2);
        alpha_exp_cl = nan; cl_exp = nan; alpha_exp_cd = nan; cd_exp = nan; alpha_exp_cc = nan; cc_exp = nan;
        alpha_mod_cl = nan; cl_mod = nan; alpha_mod_cd = nan; cd_mod = nan; alpha_mod_cc = nan; cc_mod = nan;    
    case 26
        % Test conditions
        M = 0.4;
        a_0 = 12.45*pi/180;
        a_1h = 3.14*pi/180;
        k_h = 0.116;
        authors = {'Liiva et al. (1968)'; 'Tyler & Leishman (1992)'};
        b = 6.38*0.0254/2;
        a_inf = 325;
        ah = -1/2;
        airfoil = 'Vertol-23010';
        a_1 = 0; k = 0;
        % Experimental data
        file1 = '../Other Data/CN_Vertol23010_M_04_k_0116_a0_12_a1h_3_exp.mat';
        file2 = '../Other Data/CN_Vertol23010_M_04_k_0116_a0_12_a1h_3_mod.mat';
        file3 = '../Other Data/CM_Vertol23010_M_04_k_0116_a0_12_a1h_3_exp.mat';
        file4 = '../Other Data/CM_Vertol23010_M_04_k_0116_a0_12_a1h_3_mod.mat';
        load(file1)
        load(file2)
        load(file3)
        load(file4)
        alpha_exp_cn = CN_Vertol23010_M_04_k_0116_a0_12_a1h_3_exp(:,1);
        cn_exp = CN_Vertol23010_M_04_k_0116_a0_12_a1h_3_exp(:,2);
        alpha_exp_cm = CM_Vertol23010_M_04_k_0116_a0_12_a1h_3_exp(:,1);
        cm_exp = CM_Vertol23010_M_04_k_0116_a0_12_a1h_3_exp(:,2);
        alpha_mod_cn = CN_Vertol23010_M_04_k_0116_a0_12_a1h_3_mod(:,1);
        cn_mod = CN_Vertol23010_M_04_k_0116_a0_12_a1h_3_mod(:,2);
        alpha_mod_cm = CM_Vertol23010_M_04_k_0116_a0_12_a1h_3_mod(:,1);
        cm_mod = CM_Vertol23010_M_04_k_0116_a0_12_a1h_3_mod(:,2);
        alpha_exp_cl = nan; cl_exp = nan; alpha_exp_cd = nan; cd_exp = nan; alpha_exp_cc = nan; cc_exp = nan;
        alpha_mod_cl = nan; cl_mod = nan; alpha_mod_cd = nan; cd_mod = nan; alpha_mod_cc = nan; cc_mod = nan;    
    case 27
        % Test conditions
        M = 0.4;
        a_0 = 14.88*pi/180;
        a_1h = 3.41*pi/180;
        k_h = 0.126;
        authors = {'Liiva et al. (1968)'; 'Tyler & Leishman (1992)'};
        b = 6.38*0.0254/2;
        a_inf = 325;
        ah = -1/2;
        airfoil = 'Vertol-23010';
        a_1 = 0; k = 0;
        % Experimental data
        file1 = '../Other Data/CN_Vertol23010_M_04_k_0126_a0_15_a1h_3_exp.mat';
        file2 = '../Other Data/CN_Vertol23010_M_04_k_0126_a0_15_a1h_3_mod.mat';
        file3 = '../Other Data/CM_Vertol23010_M_04_k_0126_a0_15_a1h_3_exp.mat';
        file4 = '../Other Data/CM_Vertol23010_M_04_k_0126_a0_15_a1h_3_mod.mat';
        load(file1)
        load(file2)
        load(file3)
        load(file4)
        alpha_exp_cn = CN_Vertol23010_M_04_k_0126_a0_15_a1h_3_exp(:,1);
        cn_exp = CN_Vertol23010_M_04_k_0126_a0_15_a1h_3_exp(:,2);
        alpha_exp_cm = CM_Vertol23010_M_04_k_0126_a0_15_a1h_3_exp(:,1);
        cm_exp = CM_Vertol23010_M_04_k_0126_a0_15_a1h_3_exp(:,2);
        alpha_mod_cn = CN_Vertol23010_M_04_k_0126_a0_15_a1h_3_mod(:,1);
        cn_mod = CN_Vertol23010_M_04_k_0126_a0_15_a1h_3_mod(:,2);
        alpha_mod_cm = CM_Vertol23010_M_04_k_0126_a0_15_a1h_3_mod(:,1);
        cm_mod = CM_Vertol23010_M_04_k_0126_a0_15_a1h_3_mod(:,2);
        alpha_exp_cl = nan; cl_exp = nan; alpha_exp_cd = nan; cd_exp = nan; alpha_exp_cc = nan; cc_exp = nan;
        alpha_mod_cl = nan; cl_mod = nan; alpha_mod_cd = nan; cd_mod = nan; alpha_mod_cc = nan; cc_mod = nan;  
    case 28
        % Test conditions
        M = 0.4;
        a_0 = 12.25*pi/180;
        a_1h = 1.17*pi/180;
        k_h = 0.0675;
        authors = {'Liiva et al. (1968)'};
        b = 6.38*0.0254/2;
        a_inf = 325;
        ah = -1/2;
        airfoil = 'Vertol-23010';
        a_1 = 0; k = 0;
        % Experimental data
        file1 = '../Other Data/CN_Vertol23010_M_04_k_0068_a0_125_a1h_25_exp.mat';
        file2 = '../Other Data/CM_Vertol23010_M_04_k_0068_a0_125_a1h_25_exp.mat';
        load(file1)
        load(file2)
        alpha_mod_cn = nan; cn_mod = nan; alpha_mod_cm = nan; cm_mod = nan;
        alpha_exp_cl = nan; cl_exp = nan; alpha_exp_cd = nan; cd_exp = nan; alpha_exp_cc = nan; cc_exp = nan;
        alpha_mod_cl = nan; cl_mod = nan; alpha_mod_cd = nan; cd_mod = nan; alpha_mod_cc = nan; cc_mod = nan; 
    case 29
        % Test conditions
        M = 0.4;
        a_0 = 12.46*pi/180;
        a_1h = 2.10*pi/180;
        k_h = 0.121;
        authors = {'Liiva et al. (1968)'};
        b = 6.38*0.0254/2;
        a_inf = 325;
        ah = -1/2;
        airfoil = 'Vertol-23010';
        a_1 = 0; k = 0;
        % Experimental data
        file1 = '../Other Data/CN_Vertol23010_M_04_k_0121_a0_125_a1h_25_exp.mat';
        file2 = '../Other Data/CM_Vertol23010_M_04_k_0121_a0_125_a1h_25_exp.mat';
        load(file1)
        load(file2)
        alpha_mod_cn = nan; cn_mod = nan; alpha_mod_cm = nan; cm_mod = nan;
        alpha_exp_cl = nan; cl_exp = nan; alpha_exp_cd = nan; cd_exp = nan; alpha_exp_cc = nan; cc_exp = nan;
        alpha_mod_cl = nan; cl_mod = nan; alpha_mod_cd = nan; cd_mod = nan; alpha_mod_cc = nan; cc_mod = nan;
    case 30
        % Test conditions
        M = 0.4;
        a_0 = 14.88*pi/180;
        a_1h = 2.37*pi/180;
        k_h = 0.129;
        authors = {'Liiva et al. (1968)'};
        b = 6.38*0.0254/2;
        a_inf = 325;
        ah = -1/2;
        airfoil = 'Vertol-23010';
        a_1 = 0; k = 0;
        % Experimental data
        file1 = '../Other Data/CN_Vertol23010_M_04_k_0128_a0_15_a1h_25_exp.mat';
        file2 = '../Other Data/CM_Vertol23010_M_04_k_0128_a0_15_a1h_25_exp.mat';
        load(file1)
        load(file2)
        alpha_mod_cn = nan; cn_mod = nan; alpha_mod_cm = nan; cm_mod = nan;
        alpha_exp_cl = nan; cl_exp = nan; alpha_exp_cd = nan; cd_exp = nan; alpha_exp_cc = nan; cc_exp = nan;
        alpha_mod_cl = nan; cl_mod = nan; alpha_mod_cd = nan; cd_mod = nan; alpha_mod_cc = nan; cc_mod = nan; 
    case 31
        % Test conditions
        M = 0.2;
        a_0 = 12.36*pi/180;
        a_1h = 4.64*pi/180;
        k_h = 0.254;
        authors = {'Liiva et al. (1968)'};
        b = 6.38*0.0254/2;
        a_inf = 325;
        ah = -1/2;
        airfoil = 'NACA0012';
        a_1 = 0; k = 0;
        % Experimental data
        file1 = '../Other Data/CN_NACA0012_M_02_k_024_a0_125_plunge_exp.mat';
        file2 = '../Other Data/CM_NACA0012_M_02_k_024_a0_125_plunge_exp.mat';
        load(file1)
        load(file2)
        alpha_mod_cn = nan; cn_mod = nan; alpha_mod_cm = nan; cm_mod = nan;
        alpha_exp_cl = nan; cl_exp = nan; alpha_exp_cd = nan; cd_exp = nan; alpha_exp_cc = nan; cc_exp = nan;
        alpha_mod_cl = nan; cl_mod = nan; alpha_mod_cd = nan; cd_mod = nan; alpha_mod_cc = nan; cc_mod = nan; 
    case 32
        % Test conditions
        M = 0.2;
        a_0 = 14.67*pi/180;
        a_1h = 4.60*pi/180;
        k_h = 0.249;
        authors = {'Liiva et al. (1968)'};
        b = 6.38*0.0254/2;
        a_inf = 325;
        ah = -1/2;
        airfoil = 'NACA0012';
        a_1 = 0; k = 0;
        % Experimental data
        file1 = '../Other Data/CN_NACA0012_M_02_k_024_a0_15_plunge_exp.mat';
        file2 = '../Other Data/CM_NACA0012_M_02_k_024_a0_15_plunge_exp.mat';
        load(file1)
        load(file2)
        alpha_mod_cn = nan; cn_mod = nan; alpha_mod_cm = nan; cm_mod = nan;
        alpha_exp_cl = nan; cl_exp = nan; alpha_exp_cd = nan; cd_exp = nan; alpha_exp_cc = nan; cc_exp = nan;
        alpha_mod_cl = nan; cl_mod = nan; alpha_mod_cd = nan; cd_mod = nan; alpha_mod_cc = nan; cc_mod = nan; 
    case 33
        % NREL/OSU
        [params,data] = load_OSU(374);
        [M,U,b,a_inf,beta,a_0,a_1,k,airfoil,ah,a_1h,k_h] = struct2vars(params);
        [authors,alpha_exp_cl,alpha_exp_cm,alpha_exp_cd,alpha_exp_cn,alpha_exp_cc,cl_exp,cm_exp,cd_exp,cn_exp,cc_exp] = struct2vars(data);
        % Mohamed & Wood (2020)
        file1 = '../Other Data/S809_Mohamed_N1_cn.mat';
        file2 = '../Other Data/S809_Mohamed_N1_cm.mat';
        file3 = '../Other Data/S809_Mohamed_N1_cc.mat';
        load(file1)
        load(file2)
        load(file3)
        alpha_mod_cn = S809_Mohamed_N1_cn(:,1);
        cn_mod       = S809_Mohamed_N1_cn(:,2);
        alpha_mod_cm = S809_Mohamed_N1_cm(:,1);
        cm_mod       = S809_Mohamed_N1_cm(:,2);
        alpha_mod_cc = S809_Mohamed_N1_cc(:,1);
        cc_mod       = S809_Mohamed_N1_cc(:,2);
        alpha_mod_cl = nan; cl_mod = nan; alpha_mod_cd = nan; cd_mod = nan;
        authors{2} = 'Mohamed and Wood (2020)';
    case 34
        % NREL/OSU
        [params,data] = load_OSU(398);
        [M,U,b,a_inf,beta,a_0,a_1,k,airfoil,ah,a_1h,k_h] = struct2vars(params);
        [authors,alpha_exp_cl,alpha_exp_cm,alpha_exp_cd,alpha_exp_cn,alpha_exp_cc,cl_exp,cm_exp,cd_exp,cn_exp,cc_exp] = struct2vars(data);
        % Mohamed & Wood (2020)
        file1 = '../Other Data/S809_Mohamed_N3_cn.mat';
        file2 = '../Other Data/S809_Mohamed_N3_cm.mat';
        file3 = '../Other Data/S809_Mohamed_N3_cc.mat';
        load(file1)
        load(file2)
        load(file3)
        alpha_mod_cn = S809_Mohamed_N3_cn(:,1);
        cn_mod       = S809_Mohamed_N3_cn(:,2);
        alpha_mod_cm = S809_Mohamed_N3_cm(:,1);
        cm_mod       = S809_Mohamed_N3_cm(:,2);
        alpha_mod_cc = S809_Mohamed_N3_cc(:,1);
        cc_mod       = S809_Mohamed_N3_cc(:,2);
        alpha_mod_cl = nan; cl_mod = nan; alpha_mod_cd = nan; cd_mod = nan;
        authors{2} = 'Mohamed and Wood (2020)'; 
    case 35
        % NREL/OSU
        [params,data] = load_OSU(389);
        [M,U,b,a_inf,beta,a_0,a_1,k,airfoil,ah,a_1h,k_h] = struct2vars(params);
        [authors,alpha_exp_cl,alpha_exp_cm,alpha_exp_cd,alpha_exp_cn,alpha_exp_cc,cl_exp,cm_exp,cd_exp,cn_exp,cc_exp] = struct2vars(data);
        % Mohamed & Wood (2020)
        file1 = '../Other Data/S809_Mohamed_N8_cn.mat';
        file2 = '../Other Data/S809_Mohamed_N8_cm.mat';
        file3 = '../Other Data/S809_Mohamed_N8_cc.mat';
        load(file1)
        load(file2)
        load(file3)
        alpha_mod_cn = S809_Mohamed_N8_cn(:,1);
        cn_mod       = S809_Mohamed_N8_cn(:,2);
        alpha_mod_cm = S809_Mohamed_N8_cm(:,1);
        cm_mod       = S809_Mohamed_N8_cm(:,2);
        alpha_mod_cc = S809_Mohamed_N8_cc(:,1);
        cc_mod       = S809_Mohamed_N8_cc(:,2);
        alpha_mod_cl = nan; cl_mod = nan; alpha_mod_cd = nan; cd_mod = nan;
        authors{2} = 'Mohamed and Wood (2020)'; 
    case 36
        % NREL/OSU
        [params,data] = load_OSU(392);
        [M,U,b,a_inf,beta,a_0,a_1,k,airfoil,ah,a_1h,k_h] = struct2vars(params);
        [authors,alpha_exp_cl,alpha_exp_cm,alpha_exp_cd,alpha_exp_cn,alpha_exp_cc,cl_exp,cm_exp,cd_exp,cn_exp,cc_exp] = struct2vars(data);
        % Mohamed & Wood (2020)
        file1 = '../Other Data/S809_Mohamed_N5_cn.mat';
        file2 = '../Other Data/S809_Mohamed_N5_cm.mat';
        file3 = '../Other Data/S809_Mohamed_N5_cc.mat';
        load(file1)
        load(file2)
        load(file3)
        alpha_mod_cn = S809_Mohamed_N5_cn(:,1);
        cn_mod       = S809_Mohamed_N5_cn(:,2);
        alpha_mod_cm = S809_Mohamed_N5_cm(:,1);
        cm_mod       = S809_Mohamed_N5_cm(:,2);
        alpha_mod_cc = S809_Mohamed_N5_cc(:,1);
        cc_mod       = S809_Mohamed_N5_cc(:,2);
        alpha_mod_cl = nan; cl_mod = nan; alpha_mod_cd = nan; cd_mod = nan;
        authors{2} = 'Mohamed and Wood (2020)';     
    otherwise
        error('test case not available')
end

% Additional flow variables
U = M*a_inf;
beta = sqrt(1-M^2);

%% Set all flow and test condition variables into the params struct
params = variables2struct(struct,airfoil,ah,M,U,b,a_inf,beta,a_0,a_1,k,a_1h,k_h);

%% Set on data structure
data = variables2struct(struct,authors,alpha_exp_cl,alpha_exp_cm,alpha_exp_cd,alpha_exp_cn,alpha_exp_cc,cl_exp,cm_exp,cd_exp,cn_exp,cc_exp,alpha_mod_cl,alpha_mod_cm,alpha_mod_cd,alpha_mod_cn,alpha_mod_cc,cl_mod,cm_mod,cd_mod,cn_mod,cc_mod);

end