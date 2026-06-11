function p = ebar(var,xind,yind,nf,xscale,yscale)

[x,xbar] = barfcn(var,xind,nf,xscale);
[y,ybar] = barfcn(var,yind,nf,yscale);

plot(xbar,y,'Color',[0.7 0.7 0.7]);
plot(x,ybar,'Color',[0.7 0.7 0.7]);
p = plot(x(1),y(1),'o');

return


function [v,vbar] = barfcn(var,ind,nf,scale)

% normalize error bar
if nf
    var(:,2) = var(:,2)./sqrt(var(:,3));
end

if length(ind) == 1
    v = var(ind,1);
    dv = var(ind,2);
else
    v = [10^var(ind(1),1) 10^var(ind(2),1)];
    dv = v.*var(ind,2)';
    dv = sqrt( (dv(1)/v(2))^2 + (dv(2)*v(1)/v(2)^2)^2 );
    v = v(1)/v(2);
    dv = dv/(v*log(10));
    v = log10(v);
end

v = [v v];
vbar = v + [-2*dv 2*dv];

% convert to log scale if necessary
if scale == 0
    v = 10.^v;
    vbar = 10.^vbar;
end

return