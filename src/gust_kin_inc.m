function [wg,vg,w0,F_theta0,F2_theta0,theta0] = gust_kin_inc(wg_fun,t,U,b,lambda_g)

% Check inputs
if isempty(wg_fun)
    wg = 0; vg = 0; w0 = 0; theta0 = 0; F_theta0 = 0; F2_theta0 = 0;
    return
end

% Gust-induced normal velocity component     
w0 = wg_fun(t);

% Current gust-induced angle of attack. Notice that the equations are
% somewhat different from those provided by Leishman. For the normal force
% coefficient, the agreement is perfect, but the correct expression for the
% inertial moment coefficient (F2_theta0) could not be found.
if t >= 0
    if lambda_g > 0 % Downstream convecting gust
        % Normalized position of gust front over the chord (zeta0 = 0 at the LE and zeta0 = 1 at the TE)
        zeta0 = U/b*t/(2*lambda_g);
        % theta0 = 0 when the gust is at the LE and theta0 = pi when the gust is at the TE
        theta0 = real(acos(1-2*zeta0));
        % Expression involved in the gust-induced inertial normal coefficient
        F_theta0 = sin(2*theta0)/2-theta0;
        % Expression involved in the gust-induced inertial pitching coefficient about 1/4-chord
        if lambda_g > 1, fac = lambda_g^2; else, fac = 1/lambda_g; end
        F2_theta0 = (3/8*fac*F_theta0+1/2*sin(theta0)^3)*sign(1-lambda_g);
        % Quasi-steady gust-induced normal velocity
        wg = (1-(1-theta0/pi+sin(theta0)/pi))*w0;    
    elseif lambda_g < 0 % Upstream convecting gust
        % Normalized position of gust front over the chord (zeta0 = 0 at the LE and zeta0 = 1 at the TE)
        zeta0 = 1-U/b*t/(2*abs(lambda_g));
        % theta0 = 0 when the gust is at the LE and theta0 = pi when the gust is at the TE
        theta0 = real(acos(1-2*zeta0));
        % Expression involved in the gust-induced inertial normal coefficient
        F_theta0 = theta0-sin(2*theta0)/2;
        % Expression involved in the gust-induced inertial pitching coefficient about 1/4-chord
        if abs(lambda_g) > 1, fac = lambda_g^4; else, fac = 8*abs(lambda_g)+1/abs(lambda_g); end
        F2_theta0 = 3/8*fac*F_theta0-1/2*sin(theta0)^3;
        % Quasi-steady gust-induced normal velocity
        wg = ((pi-theta0)+sin((pi-theta0)))/pi*w0;
    else % Infinite speed gust
        theta0 = pi; F_theta0 = 0; F2_theta0 = 0; wg = w0;
    end
else
    % Gust has not reached the airfoil yet
    wg = 0;
    theta0 = 0;
    F_theta0 = 0;
    F2_theta0 = 0;
end

% Gust-induced chordwise velocity component
vg = U/lambda_g-U;

end