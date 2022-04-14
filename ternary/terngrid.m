function terngrid;
% TERNGRID - Ternary grid lines.

% grid
xa = [-0.4:0.1:0.4];
ya = zeros(size(xa));

[a,b,c] = xy2tern(xa,ya);
[xb,yb] = tern2xy(b,a,c);
[xc,yc] = tern2xy(c,b,a);

plot([xa;xb],[ya;yb],'-','LineWidth',0.25,'Color',[0.8 0.8 0.8]);
plot([xb;fliplr(xc)],[yb;fliplr(yc)],'-','LineWidth',0.25,'Color',[0.8 0.8 0.8]);
plot([xc;xa],[yc;ya],'-','LineWidth',0.25,'Color',[0.8 0.8 0.8]);

return

