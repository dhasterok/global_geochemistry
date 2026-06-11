function kseries

y = [0:0.1:10];
xu = zeros(size(y));
xm = xu;
xl = xu;

ind = y <= 3.2;
xu(ind) = (56 - 48)/(3.2 - 1.6)*(y(ind) - 3.2) + 56;
xu(~ind) = (63 - 56)/(4.0 - 3.2)*(y(~ind) - 4.0) + 63;
xu(xu + y > 100) = NaN;
xu(xu < 45) = NaN;

xm = (70 - 48)/(3.0 - 1.2)*(y - 3.0) + 70;
xm(xm + y > 100) = NaN;
xm(xm < 45) = NaN;

xl = (70 - 48)/(1.3 - 0.3)*(y - 1.3) + 70;
xl(xl + y > 100) = NaN;
xl(xl < 45) = NaN;

plot(xu,y,'r-'); hold on;
plot(xm,y,'r-');
plot(xl,y,'r-');

return