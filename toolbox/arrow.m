function arrow(x,y,varargin)

c = 1;
txt = '';
style = 'single';
while nargin > c + 1;
    switch varargin{c}
        case 'style'
            style = varargin{c+1};
            if ~any(strcmp(style,{'single','double'}))
                error('Unknown arrow style.');
            end
            c = c + 2;
        case 'text'
            txt = varargin{c+1};
            c = c + 2;
        case 'color'
            C = varargin{c+1};
            c = c + 2;
        otherwise
            error('Unknown option.');
    end
end

xl = get(gca,'XLim');
yl = get(gca,'YLim');

dx = diff(xl);
dy = diff(yl);

pos = get(gca,'Position');
R = pos(3)/pos(4);

X = x/dx;
Y = y/dy;

r = 0.05*sqrt(diff(X)^2 + diff(Y)^2);
theta = atan(diff(Y)/diff(X));
phi = pi/6;

hold on;
plot(x,y,'Color',C);

if x(1) > x(2)
    xh = x(1) - dx/R*r*[cos(theta-phi); 0; cos(theta+phi)];
else
    xh = x(1) + dx/R*r*[cos(theta-phi); 0; cos(theta+phi)];
end
if y(1) > y(2)
    yh = y(1) - dy*R*r*[sin(theta-phi); 0; sin(theta+phi)];
else
    yh = y(1) + dy*R*r*[sin(theta-phi); 0; sin(theta+phi)];
end
plot(xh,yh,'Color',C);

if strcmp(style,'double')
    if x(2) > x(1)
        xh = x(2) - dx/R*r*[cos(theta-phi); 0; cos(theta+phi)];
    else
        xh = x(2) + dx/R*r*[cos(theta-phi); 0; cos(theta+phi)];
    end
    if y(2) > y(1)
        yh = y(2) - dy*R*r*[sin(theta-phi); 0; sin(theta+phi)];
    else
        yh = y(2) + dy*R*r*[sin(theta-phi); 0; sin(theta+phi)];
    end
end
plot(xh,yh,'Color',C);

if ~strcmp(txt,'')
    if strcmp(style,'single')
        text(x(2),y(2),txt);
    elseif strcmp(style,'double')
        text(x(1) + (x(2) + x(1))/2, y(1) + (y(2) + y(1))/2,txt);
    end
end

return