function h2ofig(plutonic,volcanic);

%pH2O = nan(size(plutonic.H2O_PLUS));
%ind = ~strcmp('ocean',plutonic.COUNTRY);
%pH2O(ind) = plutonic.H2O_PLUS(ind);
pH2O = log10(plutonic.H2O_PLUS);

%vH2O = nan(size(volcanic.H2O_PLUS));
%ind = ~strcmp('ocean',volcanic.COUNTRY);
%vH2O(ind) = volcanic.H2O_PLUS(ind);
vH2O = log10(volcanic.H2O_PLUS);

figure;
subplot(221);
hold on;
plot(pH2O,log10(plutonic.heat_production),'.');
ind = strcmp('granite',plutonic.rock_type);
plot(pH2O(ind),log10(plutonic.heat_production(ind)),'.');
ind = strcmp('subalkalic gabbro',plutonic.rock_type) & plutonic.heat_production > 0 & plutonic.H2O_PLUS > 0;
plot(pH2O(ind),log10(plutonic.heat_production(ind)),'.');
hpax([-2 2],'x');
xlabel('H_2O (ppm)');
%xlim([0 10]);
%xlabel('H_2O (wt.%)');
hpax([-3 2]);
legend('plutonic','granite','gabbro');
axis square;



subplot(222);
hold on;
plot(vH2O,log10(volcanic.heat_production),'.');
ind = strcmp('rhyolite',volcanic.rock_type);
plot(vH2O(ind),log10(volcanic.heat_production(ind)),'.');
ind = strcmp('subalkalic basalt',volcanic.rock_type);
plot(vH2O(ind),log10(volcanic.heat_production(ind)),'.');
hpax([-2 2],'x');
xlabel('H_2O (ppm)');
%xlim([0 10]);
%xlabel('H_2O (wt.%)');
hpax([-3 2]);
legend('volcanic','rhyolite','basalt');
axis square;

subplot(223);
hold on;
plot(log10(abs(plutonic.CE_PPM)),log10(plutonic.heat_production),'.');
ind = strcmp('granite',plutonic.rock_type);
plot(log10(abs(plutonic.CE_PPM(ind))),log10(plutonic.heat_production(ind)),'.');
ind = strcmp('subalkalic gabbro',plutonic.rock_type) & plutonic.heat_production > 0 & plutonic.CE_PPM > 0;
plot(log10(plutonic.CE_PPM(ind)),log10(plutonic.heat_production(ind)),'.');
hpax([-1 4],'x');
xlabel('Ce (ppm)');
hpax([-3 2]);
axis square;

%ind = plutonic.heat_production > 0 & plutonic.CE_PPM > 0;
%X = [-1 4];
%p = polyfit(log10(abs(plutonic.CE_PPM(ind))),log10(plutonic.heat_production(ind)),1);
%Y = polyval(p,X);
%p = robustfit(log10(abs(plutonic.CE_PPM(ind))),log10(plutonic.heat_production(ind)))
%Y = p(1) + X*p(2);
%plot(X,Y,'-');
%subplot(221); hold on;
%plot(X-1,Y,'-');

subplot(224);
hold on;
plot(log10(abs(volcanic.CE_PPM)),log10(volcanic.heat_production),'.');
ind = strcmp('rhyolite',volcanic.rock_type);
plot(log10(abs(volcanic.CE_PPM(ind))),log10(volcanic.heat_production(ind)),'.');
ind = strcmp('subalkalic basalt',volcanic.rock_type)& volcanic.heat_production > 0 & volcanic.CE_PPM > 0;
plot(log10(abs(volcanic.CE_PPM(ind))),log10(volcanic.heat_production(ind)),'.');
hpax([-1 4],'x');
xlabel('Ce (ppm)');
hpax([-3 2]);
axis square;

%ind = volcanic.heat_production > 0 & volcanic.CE_PPM > 0;
%X = [-1 4];
%p = polyfit(log10(abs(plutonic.CE_PPM(ind))),log10(plutonic.heat_production(ind)),1);
%Y = polyval(p,X);
%p = robustfit(log10(abs(volcanic.CE_PPM(ind))),log10(volcanic.heat_production(ind)))
%Y = p(1) + X*p(2);
%plot(X,Y,'-');
%subplot(222); hold on;
%plot(X-1,Y,'-');


return