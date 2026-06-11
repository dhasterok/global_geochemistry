function hpVprov_analysis(model)

model = sortrows(model,1);
plate = unique(model.plate);

% provinces with lots of data
%nind = ~cind & model.N > 1000;
gind = strcmp(model.plate,'global');
cind = strcmp(model.plate,'continental');
oind = strcmp(model.plate,'oceanic');

pind = ~(gind | cind | oind);

nind = pind & model.r_value.^2 > 0.64 & model.N > 250 & model.m0 < 0;

colors = [117 76 41; 190 30 45; 43 57 144]/255;

% r^2 versus N
figure;
c = winter(length(plate)-1);
subplot(1,3,1:2);
Nmax = max(model.N);
h(1) = plot([0 Nmax],nanmean(model.r_value(nind).^2)*[1 1],'k--');
ltxt{1} = 'well-sampled';
hold on;
    
j = 1;
ic = 1;
for i = 1:length(plate)
    ltxt{i+1} = plate{i};
    switch plate{i}
        case {'continental','oceanic','global'}
            ind = strcmp(model.plate,plate{i});
            h(i+1) = plot([0 Nmax],model.r_value(ind)^2*[1 1],'-','Color',colors(ic,:));      
            ic = ic + 1;
        otherwise
            ind = strcmp(model.plate,plate{i}) & strcmp(model.crust_type{i},'oceanic');
            s = scatter(model.N(ind),model.r_value(ind).^2,[],c(j,:));
            s.MarkerFaceColor = [1 1 1];
            
            ind = strcmp(model.plate,plate{i}) & strcmp(model.crust_type{i},'continental');
            h(i+1) = scatter(model.N(ind),model.r_value(ind).^2,[],c(j,:),'filled');
            
            j = j + 1;
            hold on;
    end
end
legend(h,ltxt,'Location','southeast');
xlabel('Number of analyses');
ylabel('r^2 value');
set(gca,'Box','on');
ylim([0 1]);
xlim([0 round(max(model.N(nind)),-3)]);

subplot(133);
yl = [0 120];
h = histogram(model.r_value(pind).^2,'DisplayStyle','stairs','BinWidth',0.025);
h.EdgeColor = colors(2,:);
hold on;
h = histogram(model.r_value(nind  & strcmp(model.crust_type,'continental')).^2,'DisplayStyle','stairs','BinWidth',0.025);
h.EdgeColor = colors(1,:);
h = histogram(model.r_value(nind).^2,'DisplayStyle','stairs','BinWidth',0.025);
h.EdgeColor = colors(3,:);
plot(model.r_value(gind)^2*[1 1],yl,'r-');
plot(mean(model.r_value(nind).^2,'omitnan')*[1 1],yl,'k--');
xlabel('Correlation Coefficient');
xlim([0 1]);
set(gca,'View',[90 -90]);
%golden;

% slope and intercept versus N
figure;
subplot(121);
fill([0 Nmax Nmax 0], ...
    model.m0(gind) + [-model.m0_a95(gind) -model.m0_a95(gind) ...
    model.m0_a95(gind) model.m0_a95(gind)], ...
    [0.8 0.8 0.8]);
hold on;
plot([0 Nmax],model.m0(gind)*[1 1],'r-');
plot([0 Nmax],mean(model.m0(nind),'omitnan')*[1 1],'k--');
for i = 1:length(plate)
    if strcmp(plate{i},'global') || strcmp(plate{i},'continental') || strcmp(plate{i},'oceanic')
        continue;
    end
    ind = strcmp(model.plate,plate{i});
    h = scatterbar(model.N(ind),model.m0(ind),model.m0_a95(ind),'y');
end
%legend(plate);
xlabel('Number of analyses');
ylabel('slope');
set(gca,'Box','on');
ylim([-15 15]);


subplot(122);
fill([0 Nmax Nmax 0], ...
    model.m1(gind) + [-model.m1_a95(gind) -model.m1_a95(gind) ...
    model.m1_a95(gind) model.m1_a95(gind)], ...
    [0.8 0.8 0.8]);
hold on;
plot([0 Nmax],model.m1(gind)*[1 1],'r-');
plot([0 Nmax],mean(model.m1(nind),'omitnan')*[1 1],'k--');
for i = 1:length(plate)
    if strcmp(plate{i},'global') || strcmp(plate{i},'continental') || strcmp(plate{i},'oceanic')
        continue;
    end
    ind = strcmp(model.plate,plate{i});
    h = scatterbar(model.N(ind),model.m1(ind),model.m1_a95(ind),'y');
    hold on;
end
%legend(plate);
xlabel('Number of analyses');
ylabel('intercept');
set(gca,'Box','on');
ylim([-5 5]);


% histograms
figure;
subplot(211);
fill(model.m0(gind) + [-model.m0_a95(gind) -model.m0_a95(gind) ...
    model.m0_a95(gind) model.m0_a95(gind)], [yl fliplr(yl)], ...
    [0.8 0.8 0.8]);
hold on;
h = histogram(model.m0(pind),'DisplayStyle','stairs','BinWidth',0.5);
h.EdgeColor = colors(2,:);
h = histogram(model.m0(nind & strcmp(model.crust_type,'continental')),'DisplayStyle','stairs','BinWidth',0.5);
h.EdgeColor = colors(1,:);
h = histogram(model.m0(nind),'DisplayStyle','stairs','BinWidth',0.5);
h.EdgeColor = colors(3,:);
plot(model.m0(gind)*[1 1],yl,'r-');
plot(mean(model.m0(nind),'omitnan')*[1 1],yl,'k--');
xlabel('Slope [log_{10}(\muW m^{-3}) s km^{-1}]');
golden;

subplot(212);
fill(model.m1(gind) + [-model.m1_a95(gind) -model.m1_a95(gind) ...
    model.m1_a95(gind) model.m1_a95(gind)], [yl fliplr(yl)], ...
    [0.8 0.8 0.8]);
hold on;
h = histogram(model.m1(pind),'DisplayStyle','stairs','BinWidth',0.2);
h.EdgeColor = colors(2,:);
h = histogram(model.m1(nind & strcmp(model.crust_type,'continental')),'DisplayStyle','stairs','BinWidth',0.2);
h.EdgeColor = colors(1,:);
h = histogram(model.m1(nind),'DisplayStyle','stairs','BinWidth',0.2);
h.EdgeColor = colors(3,:);
plot(model.m1(gind)*[1 1],yl,'r-');
plot(mean(model.m1(nind),'omitnan')*[1 1],yl,'k--');
xlabel('Intercept [log_{10}(\muW m^{-3})]');
golden;
xlim([-2 3]);

return
