function [t,wg,wgdot] = get_gust_fom_PSD(gust_profile,U,tf,options) 

% Unpack options
dt = options.dt;            % Time resolution [s]
L = options.L;              % Length scale = 2500 ft
fmax = options.fmax;        % Maximum frequency [Hz]
sigma_w = options.sigma_w;  % RMS of gust velocity [m/s]

% Define frequency vector for input PSD
omega_min = 0;                          % Minimum frequency
df = min([1/tf; 1/fmax]);               % Default resolution [Hz] - This ensures a mean of zero and a stable value on the RMS of the gust velocity profile
d_omega = df*(2*pi);                    % Frequency resolution [rad/s]
omega_max = fmax*(2*pi);                % Maximum frequency [rad/s]
omega = omega_min:d_omega:omega_max;    % Frequency vector [rad/s]

% Get PSD [units of (m/s)^2/(rad/s)]
if contains(gust_profile,"vK")
    PHI = sigma_w^2*L/pi/U*((1+8/3*(1.339*L*omega/U).^2)./(1+(1.339*L*omega/U).^2).^(11/6));
elseif contains(gust_profile,"D")
    PHI = sigma_w^2*L/pi/U*((1+3*L^2*(omega/U).^2)./(1+L^2*(omega/U).^2).^2);
end

% Time vector
t = (0:dt:tf)'; Nt = length(t);

% Random phase vector
psi = 2*pi*rand(length(omega),1);

% Gust velocity vector
tmp = sqrt(2*repmat(PHI,Nt,1)*d_omega).*cos(t*omega+repmat(psi',Nt,1));
wg = sum(tmp,2)-tmp(:,1);

% Gust velocity rate vector
wgdot = get_wgdot_PSD(t,wg,dt);

end