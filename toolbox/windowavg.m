function [xm,ym] = windowavg(x,w)

x = sortrows(x);
n = ceil(100*(max(x(:,1)) - min(x(:,1)))/w);
xm = linspace(min(x(:,1)),max(x(:,1)),n+1)';

ym = zeros(length(xm)-1,1);
for i = 1:n
    ind = xm(i)-w/2 <= x(:,1) & x(:,1) <= xm(i+1)+w/2;
    ym(i) = median(x(ind,2));
end
xm = midpt(xm);

return