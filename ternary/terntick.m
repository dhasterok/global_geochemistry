function terntick;
% TERNTICK - Ternary tick marks.

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

return
