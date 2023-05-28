function P = hemmersly(N,D,LB,UB,varargin)

% INPUTS:
% N = number of particles
% D = number of dimensions (or optimization variables)
% LB and UB = lower and upper bounds for each dimension
% OPTIONAL INPUTS:
% d1 = first dimension to be allocated
% OUTPUTS:
% P = position of the particles (N x D)

% Ref: Wong, Tien-Tsin, Wai-Shing Luk, and Pheng-Ann Heng. "Sampling with 
% Hammersley and Halton points." Journal of Graphics Tools - 1997

% Handle inputs
if ~isempty(varargin)
    d1 = varargin{1};
else
    d1 = randi([1,D]);
end

% Prime numbers list
PRIMES = primes(600); 

% Initialize output
P = nan(N,D);

% Loop over particles
for s=1:N
    % Assign first dimension
    P(s,d1) = LB(d1)+(UB(d1)-LB(d1))*(s-1)/(N-1);  
    % Loop over remaining dimensions
    dd = d1;
    for d = 2:D
        % Index for prime number
        dd = dd+1;
        if dd > D
            dd = 1; % Restart index when it goes beyond the number of dimensions
        end
        % Initialize
        p = PRIMES(dd);     % p
        pp = p;             % p_prime
        kp = s;             % k_prime
        phi = 0;            % phi function
        % Algorithm for phi
        while kp > 0
            a = mod(kp,p);
            phi = phi + a/pp;
            kp = floor(kp/p);
            pp = pp*p;
        end
        % Set dimension value for current particle
        P(s,dd) = LB(dd) + phi*(UB(dd)-LB(dd));
    end
end

end