function [m,merr] = wllsq(x,y,sigy)

ind = isnan(sigy);
sigy(ind) = 0.15*y(ind);

ind = ~isnan(y);
x = x(ind);
y = y(ind);
sigy = sigy(ind);

x = x(:);
y = y(:);

A = [ones(size(x)) x];
W = diag(1./sigy.^2);

X = inv(A'*W*A);
m = X*A'*W*y;

S = sum(sigy.^-2.*(y - A*m).^2);

merr = S/(length(x)-1).*diag(X);

return
