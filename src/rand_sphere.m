function r = rand_sphere(size,D,radii)

R = randn(size,D);                             % size x D Gaussian random number

r = R./vecnorm(R,2,2).*rand(size,1).*radii;    % Normalized Gaussian distributed points inside a D-sphere

end