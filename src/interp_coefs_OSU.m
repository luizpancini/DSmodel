function [time_interp,alpha_interp,cn_interp,cm_interp,cc_interp,cl_interp,cd_interp] = interp_coefs_OSU(data,N_interp,interp_method)

% Unpack experimental data
time_cycles = data.time_cycles;
alpha_cycles = data.alpha_cycles;
cn_cycles = data.cn_cycles;
cm_cycles = data.cm_cycles;
cc_cycles = data.cc_cycles;
cl_cycles = data.cl_cycles;
cd_cycles = data.cd_cycles;

%% Cycle angle (time) vector for interpolation
time_interp = linspace(-90,270,N_interp);

%% Get interpolated experimental values separated by cycle
% Number of complete cycles
n_cycles = length(time_cycles); 
% Initiliaze interpolation arrays for each cycle
cn_interp_cyc = nan(N_interp,n_cycles);
cm_interp_cyc = cn_interp_cyc;
cc_interp_cyc = cn_interp_cyc;
cl_interp_cyc = cn_interp_cyc;
cd_interp_cyc = cn_interp_cyc;
alpha_interp_cyc = cn_interp_cyc;
% Loop over cycles
for n=1:n_cycles
    % Interpolate coefficients
    cn_interp_cyc(:,n) = interp1(time_cycles{n},cn_cycles{n},time_interp,interp_method);
    cm_interp_cyc(:,n) = interp1(time_cycles{n},cm_cycles{n},time_interp,interp_method);
    cc_interp_cyc(:,n) = interp1(time_cycles{n},cc_cycles{n},time_interp,interp_method);
    cl_interp_cyc(:,n) = interp1(time_cycles{n},cl_cycles{n},time_interp,interp_method);
    cd_interp_cyc(:,n) = interp1(time_cycles{n},cd_cycles{n},time_interp,interp_method);
    % Interpolate AoA
    alpha_interp_cyc(:,n) = interp1(time_cycles{n},alpha_cycles{n},time_interp,interp_method);
end
% Take mean values over cycles
cn_interp = mean(cn_interp_cyc,2,'omitnan')';
cm_interp = mean(cm_interp_cyc,2,'omitnan')';
cc_interp = mean(cc_interp_cyc,2,'omitnan')';
cl_interp = mean(cl_interp_cyc,2,'omitnan')';
cd_interp = mean(cd_interp_cyc,2,'omitnan')';
alpha_interp = mean(alpha_interp_cyc,2,'omitnan')';

end