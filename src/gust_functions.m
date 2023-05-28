function [wg_fun,wgdot_fun,wgddot_fun] = gust_functions(gust_profile,options,U,b,tf)

switch gust_profile   
    case {"PSD-vK","PSD-D"}
        [wgt,wgv,wgdotv] = get_gust_fom_PSD(gust_profile,U,tf,options); 
        wg_fun = @(t) lininterp1(wgt,wgv,t);
        wgdot_fun = @(t) lininterp1(wgt,wgdotv,t);
    case "sharp-edge"
        wg_0 = options.wg_0; H = options.H;
        wg_fun = @(t) wg_0*(heaviside(t)-heaviside(t-2*H*b/U)); 
        wgdot_fun = @(t) 0;
        wgddot_fun = @(t) 0;
    case "1-cos"
        wg_0 = options.wg_0; H = options.H;
        wg_fun = @(t) 1/2*wg_0*(1-cos(2*pi*t/(2*H*b/U))).*(heaviside(t)-heaviside(t-2*H*b/U)); 
        wgdot_fun = @(t) (wg_0*(dirac((2*H*b - U*t)/U) - dirac(t))*(cos((pi*U*t)/(H*b)) - 1))/2 - (U*wg_0*pi*sin((pi*U*t)/(H*b)).*(heaviside(-(2*H*b - U*t)/U) - heaviside(t)))/(2*H*b);
        wgddot_fun = @(t) - (wg_0*(dirac(1, (2*H*b - U*t)/U) + dirac(1, t))*(cos((pi*U*t)/(H*b)) - 1))/2 - (U^2*wg_0*pi^2*cos((pi*U*t)/(H*b))*(heaviside(-(2*H*b - U*t)/U) - heaviside(t)))/(2*H^2*b^2) - (U*wg_0*pi*sin((pi*U*t)/(H*b))*(dirac((2*H*b - U*t)/U) - dirac(t)))/(H*b);
    case "sine"
        wg_0 = options.wg_0; f_sin = options.f_sin; H = options.H;
        wg_fun = @(t) wg_0*sin(f_sin*2*pi*t).*(heaviside(t)-heaviside(t-H*b/U)); 
        wgdot_fun = @(t) 2*f_sin*wg_0*pi*cos(2*pi*f_sin*t).*(heaviside(t)-heaviside(t-2*H*b/U));
        wgddot_fun = @(t) wg_0*sin(2*pi*f_sin*t)*(dirac(1, (H*b - U*t)/U) + dirac(1, t)) - 4*f_sin*wg_0*pi*cos(2*pi*f_sin*t)*(dirac((H*b - U*t)/U) - dirac(t)) + 4*f_sin^2*wg_0*pi^2*sin(2*pi*f_sin*t)*(heaviside(-(H*b - U*t)/U) - heaviside(t));
    case "swept-sine"
        wg_0 = options.wg_0; f_i = options.f_i; f_f = options.f_f; tg_end = options.tg_end;
        wg_fun = @(t) wg_0*sin((f_i+(f_f-f_i)*t/tg_end)*2*pi*t); 
        wgdot_fun = @(t) wg_0*cos(t*pi*(2*f_i + (2*t*(f_f - f_i))/tg_end))*(pi*(2*f_i + (2*t*(f_f - f_i))/tg_end) + (2*t*pi*(f_f - f_i))/tg_end);
        wgddot_fun = @(t) (4*wg_0*pi*cos(t*pi*(2*f_i + (2*t*(f_f - f_i))/tg_end))*(f_f - f_i))/tg_end - wg_0*sin(t*pi*(2*f_i + (2*t*(f_f - f_i))/tg_end))*(pi*(2*f_i + (2*t*(f_f - f_i))/tg_end) + (2*t*pi*(f_f - f_i))/tg_end)^2;
    otherwise
        wg_fun = []; wgdot_fun = []; wgddot_fun = [];
end

end