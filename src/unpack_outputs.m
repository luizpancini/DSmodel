function varargout = unpack_outputs(outputs,model)

switch model
    case {"BL","BLT"}
        t = outputs.t;
        alpha = outputs.alpha;
        alpha_plunge = outputs.alpha_plunge;
        alpha_QS = outputs.alpha_QS;
        q_QS = outputs.q_QS; 
        qR = outputs.qR;  
        R = outputs.R;  
        alpha_cr = outputs.alpha_cr; 
        theta = outputs.theta;   
        theta_min = outputs.theta_min;
        theta_max = outputs.theta_max;
        S = outputs.S;
        P = outputs.P;
        T = outputs.T;
        alpha1_n = outputs.alpha1_n;
        alpha1_m = outputs.alpha1_m;
        alpha1_c = outputs.alpha1_c;
        dalpha1_n = outputs.dalpha1_n;
        dalpha1_m = outputs.dalpha1_m;
        dalpha1_c = outputs.dalpha1_c;
        f_n = outputs.f_n;
        f_m = outputs.f_m;
        f_c = outputs.f_c;
        fprime_n = outputs.fprime_n;
        fprime_m = outputs.fprime_m;
        fprime_c = outputs.fprime_c;
        Tf_n = outputs.Tf_n;
        Tf_m = outputs.Tf_m;
        Tf_c = outputs.Tf_c;
        Ta_theta = outputs.Ta_theta;
        alpha_C = outputs.alpha_C;
        c_n = outputs.c_n;
        c_nC = outputs.c_nC;
        c_nI = outputs.c_nI;
        c_nf = outputs.c_nf;
        c_nv = outputs.c_nv;
        c_m = outputs.c_m;
        c_mC = outputs.c_mC;
        c_mI = outputs.c_mI;
        c_mf = outputs.c_mf;
        c_mv = outputs.c_mv;
        dCP = outputs.dCP;
        c_c = outputs.c_c;
        c_l = outputs.c_l;
        c_d = outputs.c_d;
        alpha_lag = outputs.alpha_lag;
        f2prime_n = outputs.f2prime_n;
        f2prime_m = outputs.f2prime_m;
        f2prime_c = outputs.f2prime_c;
        RD = outputs.RD;
        RD_theta = outputs.RD_theta;
        RD_tv0 = outputs.RD_tv0;
        f_diff_tv0 = outputs.f_diff_tv0;
        TvL_tv0 = outputs.TvL_tv0;
        f_diff_tv0_2 = outputs.f_diff_tv0_2;
        varargout = {t,alpha,alpha_plunge,alpha_QS,q_QS,qR,R,alpha_cr,theta,theta_min,theta_max,S,P,T,alpha1_n,alpha1_m,alpha1_c,dalpha1_n,dalpha1_m,dalpha1_c,f_n,f_m,f_c,fprime_n,fprime_m,fprime_c,Tf_n,Tf_m,Tf_c,Ta_theta,alpha_C,c_n,c_nC,c_nI,c_nf,c_nv,c_m,c_mC,c_mI,c_mf,c_mv,dCP,c_c,c_l,c_d,alpha_lag,f2prime_n,f2prime_m,f2prime_c,RD,RD_theta,RD_tv0,f_diff_tv0,TvL_tv0,f_diff_tv0_2};
end

end