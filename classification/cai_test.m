function cai_test(data)

ind = rockgroup(data,'all igneous') & data.heat_production > 0;

type = {'calcic','calc-alkalic','alkali-calcic','alkalic'};

figure;
subplot(121);
for i = 1:length(type);
    indt = ind & strcmp(data.frost_class2,type{i});
    plot(data.CAI(indt),log10(data.heat_production(indt)),'.');
    hold on;
end
yl = get(gca,'YLim');
plot([-2.1 -2.1],yl,'k-');
plot([0 0],yl,'k-');
plot([2.6 2.6],yl,'k-');
legend(type)
axis square;
xlabel('calcic-alkalic index (CAI)');

subplot(122);
sio2 = repmat([50:80]',[1,3]);
p = [-45.36 1.0043 -0.00427;
    -44.72 1.094 -0.00527;
    -41.86 1.12 -0.00572];
for i = 3:-1:1;
    m(:,i) = p(i,1) + p(i,2)*sio2(:,i) + p(i,3)*sio2(:,i).^2;
end
plot(sio2,m,'-');
hold on;
plot(sio2,m(:,2)-2.1,'k--',sio2,m(:,2)+2.6,'k--');
xlabel('SiO_2 (wt.%)');
ylabel('modified alkali-lime index (MALI)');
axis square;

cipw = cipwnorm(data);
fsp = {'Orthoclase','Albite','Anorthite', ...
    'Nepheline', ...%'Kaliophilite','Leucite', ...
    'Quartz'};

figure;
total = sum(cipw{:,fsp},2);
for i = 1:length(type)
    indt = ind & strcmp(data.frost_class2,type{i});
    subplot(length(type),1,i); hold on;
    for j = 1:length(fsp)
        x = cipw{indt,fsp{j}};%./total(indt);
        x(x == 0) = NaN;
        histogram(x,'BinWidth',1,'DisplayStyle','stairs');
    end
    title(type{i});
    legend(fsp);
    set(gca,'Box','on');
end
xlabel('percentage of all normative feldspa(ars/thoids)');


edges = [-inf -2.1 0 2.6 inf];

figure;
total = sum(cipw{:,fsp},2);
for i = 1:length(type)
    indt = ind & edges(i) < data.CAI & data.CAI <= edges(i+1);
    subplot(length(type),1,i); hold on;
    for j = 1:length(fsp)
        x = cipw{indt,fsp{j}};%./total(indt);
        x(x == 0) = NaN;
        histogram(x,'BinWidth',1,'DisplayStyle','stairs');
    end
    title([type{i},', CAI (',num2str(edges(i)),',',num2str(edges(i+1)),']']);
    legend(fsp);
    set(gca,'Box','on');
end
xlabel('percentage of all normative feldspa(ars/thoids)');

return