function arrow(x,y,varargin)

xl = get(gca,'XLim');
yl = get(gca,'YLim');

dx = diff(xl);
dy = diff(yl);

pos = get(gca,'Position');
R = pos(3)/pos(4);

X = x/dx;
Y = y/dy;

r = 0.1*sqrt(diff(X)^2 + diff(Y)^2);
theta = atan(diff(Y)/diff(X));
phi = pi/4;

xh = x(1) + dx/R*r*[cos(theta-phi); 0; cos(theta+phi)];
yh = y(1) + dy*R*r*[sin(theta-phi); 0; sin(theta+phi)];

hold on;
plot(x,y,'r-');
plot(xh,yh,'r-');


if nargin == 3
    text(x(2),y(2),varargin{1});
end

return