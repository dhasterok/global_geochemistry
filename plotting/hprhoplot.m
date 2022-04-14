function [m,merr,C] = hprhoplot(rho,A,Rho,colour);

for i = 1:length(Rho)-1
    ind = find(Rho(i) <= rho & rho < Rho(i+1));

    Atemp{i} = A(ind);
    rhotemp{i} = rho(ind);
end

N = 0;
for i = 1:length(Atemp)
    N = N + length(Atemp{i});
end

[X,Y] = whisker(rhotemp,Atemp,'Color',colour,'Scale','log');

hold on;
ind = ~isnan(Y(:,3)) & 2675 <= X(:,3) & X(:,3) <= 3050;
tmp = corrcoef(X(ind,3),Y(ind,3));
C = tmp(1,2);

[m,s] = polyfit(X(ind,3)-2600,Y(ind,3),1);

invR = inv(s.R);
merr = sqrt(diag(invR*invR')*s.normr^2/s.df);

p = plot(Rho,polyval(m,Rho-2600),'-');
set(p,'Color',colour);

hpax([-2 1],'y');
xlabel('Estimated Density [kg m^{-3}]');
xlim([2500 3400]);
set(gca,'Box','on');
golden;

return

