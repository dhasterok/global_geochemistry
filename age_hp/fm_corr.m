function fm_corr(data,agediv,sio2split);

Q = [0.05 0.25 0.50 0.75 0.95];
errfloor = 0.05;
agemid = midpt(agediv);

for i = 1:length(agediv)-1
    ind = agediv(i) <= data.avg_age & data.avg_age < agediv(i+1);

    indf = ind & data.sio2 <= sio2split;
    mafic.n(i,1) = sum(indf);
    mafic.age(i,:) = quantile(data.avg_age(indf),Q);
    mafic.hp(i,:) = quantile(log10(data.heat_production(indf)),Q);
    
    indf = ind & data.sio2 > sio2split;
    felsic.n(i,1) = sum(indf);
    felsic.age(i,:) = quantile(data.avg_age(indf),Q);
    felsic.hp(i,:) = quantile(log10(data.heat_production(indf)),Q);
end
mafic.errbar(:,1) = mafic.hp(:,3) + ...
    (mafic.hp(:,5) - mafic.hp(:,3))./sqrt(mafic.n);
mafic.errbar(:,2) = mafic.hp(:,3) - ...
    (mafic.hp(:,3) - mafic.hp(:,1))./sqrt(mafic.n);
felsic.errbar(:,1) = felsic.hp(:,3) + ...
    (felsic.hp(:,5) - felsic.hp(:,3))./sqrt(felsic.n);
felsic.errbar(:,2) = felsic.hp(:,3) - ...
    (felsic.hp(:,3) - felsic.hp(:,1))./sqrt(felsic.n);

hold on;
ind = agemid < 3000;
plot([-2 1],[-2 1]+0.7,'-');
[m,r,rsq] = wllsq(mafic.hp(ind,3),felsic.hp(ind,3), ...
    max(diff(mafic.errbar(ind,:),1,2),errfloor*ones([sum(ind),1])), ...
    max(diff(felsic.errbar(ind,:),1,2),errfloor*ones([sum(ind),1])));
plot([-2 1],m(1,1)*[-2 1] + m(2,1),'-');

p = polyfit(mafic.hp(3,:), felsic.hp(3,:),1);
plot([-2 1],polyval(p,[-2 1]),'-');
str = sprintf('A(t) = %f t + %f',p(1),p(2));
t = text(-1.2,0.6,str);

p = plot(mafic.errbar',repmat(felsic.hp(:,3),[1 2])','-');
set(p,'Color',[0.6 0.6 0.6]);
p = plot(repmat(mafic.hp(:,3),[1 2])',felsic.errbar','-');
set(p,'Color',[0.6 0.6 0.6]);

scatter(mafic.hp(:,3),felsic.hp(:,3),24,agemid,'filled');

[map,cb] = gtspal('contrast.xlsx');

for i = 1:length(agemid)
    text(mafic.hp(i,3)+0.02,felsic.hp(i,3)+0.02,num2str(agemid(i)));
end

axis square;
set(gca,'Box','on');
hpax([-2 1],'x');
xlabel('Mafic heat production [\muW m^{-3}]');
hpax([-1 1],'y');
ylabel('Felsic heat production [\muW m^{-3}]');
axis(log10([0.01 3 0.3 10]));
r = corrcoef(mafic.hp(3,:), felsic.hp(3,:));
t = text(-1.2,0.8,['C = ',num2str(r(1,2))]);

return