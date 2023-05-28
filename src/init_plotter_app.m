function [axes_size,lw,ms,xlim_vec] = init_plotter_app(a_0,a_1)

set(0,'DefaultTextInterpreter','tex')
set(0,'DefaultLegendInterpreter','tex')
axes_size = 16;
lw = 1;
ms = 5;

% Determine axes limits
if a_0+a_1 > a_0-a_1
    xlim_vec = [(a_0-a_1)*180/pi (a_0+a_1)*180/pi];
else
    xlim_vec = [(a_0+a_1)*180/pi (a_0-a_1)*180/pi];
end