function agebin = hpvsage(data,age_div,div_value)
% hpvsage(data,age_div,sio2_div)
% Plots the raw age vs. hp data and associated subplots
% Division between 'mafic' and 'felsic' wt% sio2: div_value

% Div up the data indices in a structure
for i = 1:length(age_div)-1
    ind = find(age_div(i) <= data.avg_age & data.avg_age < age_div(i+1));
    agebin.n(i) = length(ind);
    agebin.ind{i} = ind;
    avg_age{i} = data.avg_age(ind);
    heat_production{i} = data.heat_production(ind);
    field_val{i} = data{ind,'sio2'};
end

% Plot data
figure()
[agebin.Qage,agebin.Qhp] = whisker(avg_age,heat_production,'color',[0.5 0.5 0.5],'scale','linear');
%[agebin.Qage,agebin.Qhp] = whisker(avg_age,heat_production,'scale','linear');

% Raw data - standard scale
%--------------------------
subplot(8,1,1:2);
sqwavefill(agebin.Qhp,agebin.Qage(:,3),age_div,[0,0,0])
xlim([age_div(1) age_div(end)]);
set(gca,'XTickLabel',{});
set(gca,'Box','on');
ylim([0 ceil(max(agebin.Qhp(:,5)))]);
set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);
ylabel('A [\muW m^{-3}]','FontSize',10)
set(gca,'XAxisLocation','top');
xlabel('Age [Ga]');
set(gca,'XTick',[0:500:4000],'XTickLabel',[0:0.5:4]);


% Raw data - log scale
%--------------------------
subplot(8,1,3:4);
[agebin.Qage,agebin.Qhp] = whisker(avg_age,heat_production,'color',[0.5 0.5 0.5],'scale','log');
%[agebin.Qage,agebin.Qhp] = whisker(avg_age,heat_production,'scale','log');
plot(0,0)
sqwavefill(agebin.Qhp,agebin.Qage(:,3),age_div,[0,0,0])
xlim([age_div(1) age_div(end)]);
set(gca,'Box','on');
set(gca,'XTickLabel',{});
hpax([floor(min(agebin.Qhp(:,1))) ceil(max(agebin.Qhp(:,5)))]);
set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);
ylabel('A [\muW m^{-3}]','FontSize',10)

[maxval,maxind]=max(agebin.Qhp(:,3));
[minval,minind]=min(agebin.Qhp(:,3));

fprintf('Max median: %f at %f\n',10^maxval,agebin.Qage(maxind,3))
fprintf('Min median: %f at %f\n',10^minval,agebin.Qage(minind,3))


subplot(8,1,5);
histogram(data.avg_age,'BinEdges',age_div);
ylabel('N data');


% Cloud City SiO2 plots
%-----------------------
subplot(8,1,6); hold on;
for i = 1:length(field_val)
    
    if isempty(field_val{i})
        continue;
    end
    
    N = (age_div(i+1) - age_div(i))/2;
    %cloudcity(field_val{i},N,agebin.Qage(i,3),'y');
    cloudcity(field_val{i},'Scale',N,'Shift',(age_div(i+1) + age_div(i))/2);
end

%Plot dashed line at 60
plot(min(age_div):100:max(age_div),div_value.*ones(1,length(min(age_div):100:max(age_div))),'color','r','linestyle','--')

%Plots mafic/felsic text box
text(4100,48,'Mafic','rotation',90,'FontSize',8,'HorizontalAlignment','right')
text(4100,52,'Felsic','rotation',90,'FontSize',8,'HorizontalAlignment','left')
hold off
xlim([age_div(1) age_div(end)]);
set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8,...
        'XTickLabel',{});
ylabel('SiO_2 [wt %]','FontSize',10)
set(gca,'Box','on');


% Supercontinent cycle
%---------------------------
subplot(8,1,8);
hold on
orogen_hist(0,4000);
hold off


% Relative felsic/mafic percentage
%---------------------------
subplot(8,1,7)
for i = 1:length(age_div)-1
    felsic_num(i) = length(find(age_div(i) <= data.avg_age & data.avg_age < age_div(i+1) ...
        & data.sio2>div_value & data.sio2<=100));
    mafic_num(i) = length(find(age_div(i) <= data.avg_age & data.avg_age < age_div(i+1) ...
        & data.sio2<=div_value & data.sio2>0));
end
h = bar(age_div(1:end-1),(felsic_num-mafic_num)./(felsic_num+mafic_num),'histc');
set(h,'FaceColor',[0.7 0.7 0.7])
set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8,...
        'XTickLabel',{});
ylabel('SiO_2 bias','FontSize',10);
set(gca,'Box','on')
ylim([-1 1])


% Trend just related to relative 'felsic' or 'mafic'
% Make something a little more robust - 5wt% intervals * median sio2 hp for
% those wt% values?
% figure()
% weight_list = (felsic_num)./(felsic_num+mafic_num);
% felsic_hp = log10(nanmedian(data.heat_production(data.sio2 > div_value & data.sio2 > 0 & data.sio2 <= 100)));
% mafic_hp = log10(nanmedian(data.heat_production(data.sio2 <= div_value & data.sio2 > 0 & data.sio2 <= 100)));
% plot((age_div(2:end) + age_div(1:(end-1)))./2,weight_list.*felsic_hp + (1-weight_list).*mafic_hp,'-k')
% hpax([-2 2]);


% sio2list = [46:1:78];
% sio2hp = [];
% for i = 1:((length(sio2list))-1)
%     sio2hp(i) = nanmedian(data.heat_production(data.sio2 >= sio2list(i) & ...
%         data.sio2 < sio2list(i+1) & data.avg_age <= 1000));
% end
% figure()
% plot((sio2list(2:end)+ sio2list(1:(end-1)))./2,sio2hp,'-k')
% 
% age_hp_sio2 = [];
% for i = 1:length(age_div)-1
%     ind = age_div(i) <= data.avg_age & data.avg_age < age_div(i+1) & ...
%         data.heat_production >= 0 & data.sio2 < 85 & data.sio2 >= 45;
%     weights_sio2 = [];
%     for j = 1:(length(sio2list)-1)
%         weights_sio2(j) = length(find(data.sio2 >= sio2list(j) & ...
%             data.sio2 < sio2list(j+1) & ind))./length(find(ind));
%     end
%     age_hp_sio2(i) = sum(weights_sio2 .* sio2hp);
% end
% figure()
% plot((age_div(2:end) + age_div(1:(end-1)))./2,log10(age_hp_sio2),'-k')
% hpax([-2 2])
% 
% (age_div(2:end) + age_div(1:(end-1)))./2
% log10(age_hp_sio2)


return
