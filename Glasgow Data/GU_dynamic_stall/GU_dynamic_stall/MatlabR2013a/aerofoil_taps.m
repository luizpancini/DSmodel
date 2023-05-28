function [yp, dyp, lp] = aerofoil_taps (xs, ys, xp)
%
% Given the shape of the aerofoil with a number of coordinate pairs ('xs',
% 'ys') and the xp positions of the pressure tappings on the
% aerofoil surface, 'aerofoil_taps' gives the 'yp' coordinate of the
% tappings and the gradient of the surface 'dyp' at the pressure tappings.
% 'lp' is the length of the influence region of each pressure tapping.
%
% 'xs' and 'ys' are two columns matrices with coordinates from the leading
% edge to the trailing edge of the upper surface (first column) and the
% lower surface (second column).
%
% 'xp' from trailing edge on upper surface to trailing edge on lower
% surface. Second column is 1 if upper surface, -1 if lower surface (same
% for 'yp' and 'dyp').
% 



% % Plotting of the aerofoil coordinate pairs
% hFig = figure;
% plot(xs(:,1), ys(:,1), '+k')
% hold on
% plot(xs(:,2), ys(:,2), '+k')
% axis equal



% INTERPOLATION OF THE AEROFOIL SHAPE
% Interpolation step
step = 0.00005;

% Building the coordinate vectors from matrices (from trailing edge on
% upper side to trailing edge on lower side - anticlockwise)
xs_O = [flipud(xs(:,1)); xs(:,2)];
ys_O = [flipud(ys(:,1)); ys(:,2)];

% Checking and deleting repeated points (i.e. nose point of the profile)
if length(find((xs_O(2:end)-xs_O(1:end-1))==0)) == 1
    n_2 = find((xs_O(2:end)-xs_O(1:end-1))==0);
    xs_O = [xs_O(1:n_2); xs_O(n_2+2:end)];
    ys_O = [ys_O(1:n_2); ys_O(n_2+2:end)];
end

% The profiles coordinates are divided into 3 sections and interpolated along
% the x and y coordinates accordingly in order to achieve a better description
% Finding fo the maximum ordinate positions
max_y = find(ys_O==max(ys_O));
min_y = find(ys_O==min(ys_O));
% and arbitrary definition of the splitting positions of the profile sections
max_y = max_y+5;
min_y = min_y-5;

% Interpolation of the 3 sections of the aerofoil (for TE to max_y on upper
% surface, from max_y to min_y on lower surface and from min_y to TE)
x1 = [xs_O(1):-step:xs_O(max_y), xs_O(max_y)]';
y1 = interp1(xs_O(1:max_y), ys_O(1:max_y), x1, 'spline');

y2 = [ys_O(max_y):-step:ys_O(min_y), ys_O(min_y)]';
x2 = interp1(ys_O(max_y:min_y), xs_O(max_y:min_y), y2, 'spline');

x3 = [xs_O(min_y):step:xs_O(end), xs_O(end)]';
y3 = interp1(xs_O(min_y:end), ys_O(min_y:end), x3, 'spline');

% Checking repeated points
if x1(end)==x1(end-1)
    x1 = x1(1:end-1);
    y1 = y1(1:end-1);
end
if x2(end)==x2(end-1)
    x2 = x2(1:end-1);
    y2 = y2(1:end-1);
end
if x3(end)==x3(end-1)
    x3 = x3(1:end-1);
    y3 = y3(1:end-1);
end

% Interpolated coordiantes
x = [x1; x2(2:end-1); x3];
y = [y1; y2(2:end-1); y3];

% % Plotting the interpolated aerofoil shape
% figure(hFig)
% plot(x, y, '-k')
% % plot(x1, y1, '-r')
% % plot(x2, y2, '-b')
% % plot(x3, y3, '-r')



% CALCULATION OF THE PRESSURE TAPPING ORDINATES
yp = zeros(size(xp));
for i = 1 : length(xp)
    if xp(i,2)==1
        yp(i,1) = interp1(x(1:find(x==min(x))), y(1:find(x==min(x))), xp(i,1), 'spline');
    elseif xp(i,2)==-1
        yp(i,1) = interp1(x(find(x==min(x)):end), y(find(x==min(x)):end), xp(i,1), 'spline');
    else
        disp(i)
    end
end
yp(:,2) = xp(:,2);

% % Plotting of the pressure tappings
% figure(hFig)
% plot(xp(:,1), yp(:,1), 'or')



% CALCULATION OF THE AEROFOIL GRADIENT AT THE PRESSURE TAPPINGS
% Differentials
Dx = x(2:end)-x(1:end-1);
Dy = y(3:end)-y(1:end-2);

% Gradients of the interpolated coordinates
dydx = [(y(2)-y(1))/Dx(1); Dy./(Dx(1:end-1)+Dx(2:end)); (y(end)-y(end-1))/Dx(end)];

% Extrapolation of the gradients at the pressure tapping locations
dyp = zeros(size(xp));
for i = 1 : length(xp)
    if xp(i,2)==1 % upper surface
        dyp(i,1) = interp1(x(1:find(x==min(x))), dydx(1:find(x==min(x))), xp(i,1), 'spline');
    elseif xp(i,2)==-1 % lower surface
        dyp(i,1) = interp1(x(find(x==min(x)):end), dydx(find(x==min(x)):end), xp(i,1), 'spline');
    else
        disp(i)
    end
end
dyp(:,2) = xp(:,2);

% % Plotting of the aerofoil profile gradient
% figure(hFig)
% plot(x,dydx,'-b')
% xlim([-0.1 1.1])
% ylim([-0.6 0.6])
% % plot(xp(:,1), dyp(:,1),'ob')



% CALCULATION OF THE REGION OF INFLUENCE OF THE PRESSURE TAPPINGS

% Internal panels from coordinates mean points
% Distance along x and y of the pressure tapping mean points
lx = [ xs_O(1)-mean([xp(1,1), xp(2,1)]);...
      (xp(3:end,1)-xp(1:end-2,1)) - (xp(2:end-1,1)-xp(1:end-2,1));...
       xs_O(end)-mean([xp(end,1), xp(end-1,1)])];
ly = [ ys_O(1)-mean([yp(1,1), yp(2,1)]);...
      (yp(3:end,1)-yp(1:end-2,1)) - (yp(2:end-1,1)-yp(1:end-2,1));...
       ys_O(end)-mean([yp(end,1), yp(end-1,1)])];
% Panel coordinates
xl = [ xs_O(1); (xp(2:end,1)+xp(1:end-1,1))/2; xs_O(end) ];
yl = [ ys_O(1); (yp(2:end,1)+yp(1:end-1,1))/2; ys_O(end) ];

% plot(xl, yl, '-sr')


% Tangent panels from local gradient at the pressure tappings
% Panel coordinates
xl = zeros(length(xp)+1,1);
yl = zeros(length(xp)+1,1);
xl(1) = xs_O(1);
yl(1) = ys_O(1);
xl(end) = xs_O(end);
yl(end) = ys_O(end);
xl(2:end-1) = ( (yp(2:end,1)-yp(1:end-1,1)) + (xp(1:end-1,1).*dyp(1:end-1,1)-xp(2:end,1).*dyp(2:end,1)) ) ./ (dyp(1:end-1,1)-dyp(2:end,1));
yl(2:end-1) = yp(1:end-1,1) + dyp(1:end-1,1).*(xl(2:end-1)-xp(1:end-1,1));
% Panel lengths along x and y
lx = xl(2:end)-xl(1:end-1);
ly = yl(2:end)-yl(1:end-1);

% plot(xl, yl, '-sg')


% Panel lengths
lp = (lx.^2+ly.^2).^0.5;

% Plotting the region of influence length
% figure(hFig)
% plot(xp(:,1),lp,'-og')


