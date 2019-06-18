function [x,y] = getpoints(fig);
% GETPOINTS - Gets points when a left mouse button
% is clicked.  Stops when the right mouse button is
% pushed.

% Bring focus to this figure
figure(fig);
hold on;

c = 1;
while 1
    [xtemp,ytemp,button] = ginput(1);

    if button ~= 1
        return
    end
    plot(xtemp,ytemp,'rx');

    x(c,1) = xtemp;
    y(c,1) = ytemp;
    c = c + 1;
end

return
