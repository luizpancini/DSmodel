function gust_options = gust_profile_vars(gust_profile,gust_test_data)

switch gust_profile   
    case {"PSD-vK","PSD-D"}
        % Set default options
        gust_options.dt = 1e-4;         % Time resolution [s]
        gust_options.L = 0.3048*2500;   % Length scale = 2500 ft
        gust_options.fmax = 50;         % Maximum frequency [Hz]
        gust_options.sigma_w = 5;       % RMS of gust velocity [m/s]
        % Get any optional inputs
        if length(gust_test_data) > 1
            input_options = gust_test_data(2);
            for fn = fieldnames(input_options)'
                gust_options.(fn{1}) = input_options.(fn{1});
            end
        end
    case "sharp-edge"
        gust_options.wg_0 = gust_test_data.wg_0;
        gust_options.H = gust_test_data.H; 
    case "1-cos"
        gust_options.wg_0 = gust_test_data.wg_0;
        gust_options.H = gust_test_data.H; 
    case "sine"
        gust_options.wg_0 = gust_test_data.wg_0;
        gust_options.f_sin = gust_test_data.f_sin;
        gust_options.H = gust_test_data.H;
    case "swept-sine"
        gust_options.wg_0 = gust_test_data.wg_0;
        gust_options.f_i = gust_test_data.f_i;
        gust_options.f_f = gust_test_data.f_f;
        gust_options.tg_end = gust_test_data.tg_end;
    otherwise
        error("Gust profile '" + gust_profile + "' not available");
end

end