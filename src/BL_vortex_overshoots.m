function [c_nv,c_mv,c_cv] = BL_vortex_overshoots(tau_v,Tv_tv0,f_diff_tv0,qR_tv0,R_tv0,RD_tv0,theta_tv0,upstroke_tv0,f_diff_tv0_2,RD_tv0_2,upstroke_tv0_2,mu_v2,nu_1,nu_2,nu_3,nu_4,nu_5,chi_u,chi_d,g_v,g_v2,Tv,Tv2,Vm,Vc,Vn1,Vn2,Vn3)

%% First vortex
if tau_v <= 2*Tv 
    % Vortex strength - proportional to pitch rate
    Vs1 = Vn1*RD_tv0^nu_1*f_diff_tv0^nu_2*max([1,min([3,qR_tv0^nu_3])]); 
    % Vortex shape function
    chi_v = 1+chi_u*R_tv0*upstroke_tv0+chi_d*~upstroke_tv0;
    if tau_v/Tv <= 1/chi_v
        Vx_n = sin(pi/2*tau_v/Tv*chi_v)^(nu_4/chi_v+RD_tv0); 
    else
        stretch_fac = chi_v/(2*chi_v-1);
        Vx_n = sin(pi/2*(1+stretch_fac*(tau_v/Tv-1/chi_v)))^(nu_4*chi_v+RD_tv0);
    end
    % Airloads' coefficients
    cn_v1 = Vs1*Vx_n*sign(theta_tv0);
    cm_v1 = -Vm*cn_v1; 
else
    cn_v1 = 0;
    cm_v1 = 0;
end

% Chordwise force vortex airload
tv_c = Tv; % Time at which vortex begins to affect chordwise force
delta_Tvc = min([1/2,qR_tv0-R_tv0]); % Increment in vortex convection time for cc at high pitch rates
if tau_v > Tv && tau_v <= (2+delta_Tvc)*Tv 
    cc_v1 = -Vc*Vn1*RD_tv0^nu_1*f_diff_tv0^nu_2*min([5,qR_tv0^nu_5])*((1-cos(2*pi*(tau_v-tv_c)/((1+delta_Tvc)*Tv)))/2)*upstroke_tv0;
else
    cc_v1 = 0;
end

%% Second vortex 
tv2_i = (2+g_v)*Tv; % Time at which second vortex begins
if tau_v > tv2_i && tau_v <= tv2_i+2*Tv2 && f_diff_tv0_2 > 0.05
    % Vortex strength
    Vs2 = Vn2*RD_tv0^nu_1*f_diff_tv0^nu_2*max([1,min([3,qR_tv0^nu_3])])*(1+mu_v2*R_tv0^2*upstroke_tv0_2); % *(1+mu_v2*RD_tv0_2^4*upstroke_tv0_2) to increase strength of vortex shed while still in the upstroke at very high pitch rate
    % Vortex shape function - simple 1-cos 
    Vx_n2 = (1-cos(2*pi*(tau_v-tv2_i)/(2*Tv2)))/2; 
    % Airloads' coefficients
    cn_v2 = Vs2*Vx_n2*sign(theta_tv0);
    cm_v2 = -Vm*cn_v2*RD_tv0_2^3; % *RD_tv0_2^3 to reduce effect at low-medium pitch rates
else
    cn_v2 = 0;
    cm_v2 = 0;
end

%% Third vortex 
tv3_i = tv2_i+(2+g_v2)*Tv2; % Time at which third vortex begins
if tau_v > tv3_i && tau_v <= tv3_i+2*Tv2 && f_diff_tv0_2 > 0.1 && upstroke_tv0_2
    % Vortex strength
    Vs3 = Vn3*RD_tv0_2*RD_tv0^nu_1*f_diff_tv0^nu_2*(1+mu_v2*R_tv0^2);  
    % Vortex shape function 
    Vx_n3 = (1-cos(2*pi*(tau_v-tv3_i)/(2*Tv2)))/2;  
    % Airloads' coefficients
    cn_v3 = Vs3*Vx_n3*sign(theta_tv0); 
    cm_v3 = -Vm*cn_v3*RD_tv0_2^3; 
else
    cn_v3 = 0;
    cm_v3 = 0;
end 

%% Total vortices' contributions
c_nv = cn_v1+cn_v2+cn_v3;
c_mv = cm_v1+cm_v2+cm_v3;
c_cv = cc_v1;

end