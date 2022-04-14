function asi_vs_ai(data)

ind = data.ASI > 0 & data.ASI < inf & data.AI > 0 & data.AI < inf;

asi = [0:0.02:3];
ai = [0:0.02:3]';

%v = griddata(data.ASI(ind),data.AI(ind),log10(data.heat_production(ind)),asi,ai);

%surf(asi,ai,v);
subplot(121);
scatter(data.ASI(ind),data.AI(ind),4,log10(data.heat_production(ind)));
hold on;
plot([1.1 1.1],[0 3],'k--');
plot([1 1],[0 3],'k-');
plot([0 3],[1 1],'k-');

axis ([0.5 2 0.5 2]);
axis square;
axis xy;

xlabel('ASI');
ylabel('AI');
caxis([-0.5 1]);
colorbar;
set(gca,'Box','on');

subplot(122);
scatter(data.ASI(ind),data.AI(ind),4);%,data.avg_age(ind)/1e3);
hold on;
plot([1.1 1.1],[0 3],'k--');
plot([1 1],[0 3],'k-');
plot([0 3],[1 1],'k-');

axis ([0.5 2 0.5 2]);
axis square;
axis xy;

xlabel('ASI');
ylabel('AI');
caxis([0 3.8]);
colorbar;
set(gca,'Box','on');

figure;
type = {'calcic','calc-alkalic','alkali-calcic','alkalic'};
for i = 1:length(type)
    ind_t = ind & strcmp(data.frost_class2,type{i});
    
    plot(data.ASI(ind_t),data.AI(ind_t),'.');
    hold on;
end
plot([1.1 1.1],[0 3],'k--');
plot([1 1],[0 3],'k-');
plot([0 3],[1 1],'k-');

axis ([0.5 2 0.5 2]);
axis square;
axis xy;

xlabel('ASI');
ylabel('AI');
legend(type);
set(gca,'Box','on');


% ternary diagrams
cipw = cipwnorm(data(ind,:));

figure;
subplot(121);
hold on;
ternary('An','Ab','Or');
ternplot(cipw.Anorthite,cipw.Albite,cipw.Orthoclase,'o', ...
    {4*ones([sum(ind) 1]),log10(data.heat_production(ind))});
caxis([-0.5 1]);
colorbar;

subplot(122);
hold on;
ternary('Q','Ab+An','Or');
ternplot(cipw.Quartz,cipw.Albite + cipw.Anorthite,cipw.Orthoclase,'o', ...
    {4*ones([sum(ind) 1]),log10(data.heat_production(ind))});
caxis([-0.5 1]);
colorbar;

figure;
type = {'calcic','calc-alkalic','alkali-calcic','alkalic'};
subplot(121);
ternary('An','Ab','Or');

subplot(122);
ternary('Q','Ab+An','Or');
for i = 1:length(type)
    ind_t = ind & strcmp(data.frost_class2,type{i});
    
    cipw = cipwnorm(data(ind_t,:));
    
    subplot(121); hold on;
    t(i) = ternplot(cipw.Anorthite,cipw.Albite,cipw.Orthoclase,'.');
    
    subplot(122); hold on;
    ternplot(cipw.Quartz,cipw.Albite + cipw.Anorthite,cipw.Orthoclase,'.');
end
legend(t,type);

return