function ternaxes;
% TERNAXES - Creates ternary axes.

% half width
w = 0.5;
% vertical scale
h = 0.5/tan(pi/6);

hold on;
axis equal;

%create axes
fill([-w 0 w -w],[0 h 0 0],[1 1 1],'LineWidth',1);
set(gca,'Visible','off');

return
