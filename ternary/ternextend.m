function ternextend(z);

% half width
w = z;
u = w + 0.5*z;
% vertical scale
h = -0.5*z/tan(pi/6);

hold on;
axis equal;

%create axes
plot([-w -u u w],[0 h h 0],'k-','LineWidth',1);
set(gca,'Visible','off');


xt = [-0.4:0.1:0.4];
yt = zeros(size(xt));

xb = [-u+0.1:0.1:u-0.1];
yb = h*ones(size(xb));

xr = fliplr(w + 0.5*[0:0.1:z-0.1]);
xl = -xr;
ys = fliplr([0:0.2:0.8]*h);


plot([xl;xr],[ys;ys],'-','LineWidth',0.25,'Color',[0.8 0.8 0.8]);
plot([xl,xt;xb],[ys,yt;yb],'-','LineWidth',0.25,'Color',[0.8 0.8 0.8]);
plot([xr,fliplr(xt);fliplr(xb)],[ys,fliplr(yt);fliplr(yb)],'-','LineWidth',0.25,'Color',[0.8 0.8 0.8]);

axis([-1.346255486272525 1.346255486272525 -0.480230082488086 1.346255486272525]);

return