function [X,Y] = mavg2(x,y,w,n);

% sort by x values
z = [x,y];
z = sortrows(z);

x = z(:,1);
y = z(:,2);

X = linspace(x(1),x(end),n)';

for i = 1:n
    ind = X(i) - w < x & x< X(i) + w;
    Y(i) = mean(y(ind));
end

return