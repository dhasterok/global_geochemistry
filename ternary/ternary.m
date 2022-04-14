function ternary(varargin)
% TERNARY - creates ternary axes.

ternaxes(nargin);
terngrid(nargin);
terntick(nargin);

if nargin == 3
    axis([-1.346255486272525 1.346255486272525 -0.480230082488086 1.346255486272525]);
elseif nargin == 4
    axis([-1.346255486272525 1.346255486272525 -1.346255486272525 1.346255486272525]);
end

switch nargin
    case 0
    case 3
        ternlabel(varargin{1},varargin{2},varargin{3});
    case 4
        ternlabel(varargin{1},varargin{2},varargin{3},varargin{4});
    otherwise
        error('Incorrect number of input arguments.');
end
        

return

function ternaxes(n);

% half width
w = 0.5;
% vertical scale
h = 0.5/tan(pi/6);

hold on;
axis equal;

%create axes
plot([-w 0 w -w],[0 h 0 0],'k-','LineWidth',1);
if n == 4
    plot([-w 0 w -w],-[0 h 0 0],'k-','LineWidth',1);
end
set(gca,'Visible','off');

return


function terntick(n);

% tick marks
ya = zeros([3,9]);
ya([1 3],:) = 0.015;
xa(2,:) = [-0.4:0.1:0.4];
xa(1,:) = xa(2,:) + ya(1,:)/tan(pi/3);
xa(3,:) = xa(2,:) - ya(1,:)/tan(pi/3);

plot(xa,ya,'k-','LineWidth',1);
[a,b,c] = xy2tern(xa,ya);
[xb,yb] = tern2xy(b,a,c);
plot(xb,yb,'k-','LineWidth',1);
[xc,yc] = tern2xy(c,b,a);
plot(xc,yc,'k-','LineWidth',1);

if n == 4
    plot(xa,-ya,'k-','LineWidth',1);
    plot(xb,-yb,'k-','LineWidth',1);
    plot(xc,-yc,'k-','LineWidth',1);
end

return


function terngrid(n);
% grid
xa = [-0.4:0.1:0.4];
ya = zeros(size(xa));

[a,b,c] = xy2tern(xa,ya);
[xb,yb] = tern2xy(b,a,c);
[xc,yc] = tern2xy(c,b,a);

plot([xa;xb],[ya;yb],'-','LineWidth',0.25,'Color',[0.8 0.8 0.8]);
plot([xb;fliplr(xc)],[yb;fliplr(yc)],'-','LineWidth',0.25,'Color',[0.8 0.8 0.8]);
plot([xc;xa],[yc;ya],'-','LineWidth',0.25,'Color',[0.8 0.8 0.8]);

if n == 4
    plot([xa;xb],-[ya;yb],'-','LineWidth',0.25,'Color',[0.8 0.8 0.8]);
    plot([xb;fliplr(xc)],-[yb;fliplr(yc)],'-','LineWidth',0.25,'Color',[0.8 0.8 0.8]);
    plot([xc;xa],-[yc;ya],'-','LineWidth',0.25,'Color',[0.8 0.8 0.8]);

end

return


function ternlabel(atxt,btxt,ctxt,varargin);

w = 0.5;
h = 0.5/tan(pi/6);
d = 0.02;

text(0,h+d,atxt,'FontSize',12,'HorizontalAlignment','center','FontWeight','bold','VerticalAlignment','bottom');
text(-w-d,0,btxt,'FontSize',12,'HorizontalAlignment','right','FontWeight','bold');
text(w+d,0,ctxt,'FontSize',12,'HorizontalAlignment','left','FontWeight','bold');

if nargin == 4
    dtxt = varargin{1};
    text(0,-(h+d),dtxt,'FontSize',12,'HorizontalAlignment','center','FontWeight','bold','VerticalAlignment','top');
end

return