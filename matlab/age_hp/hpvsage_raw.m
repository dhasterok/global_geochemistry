function [agebin,mafic_ind,felsic_ind] = hpvsage_raw(data,age_div,folder)
%--------------------------------------------------------------------------
%
%                 Age vs. HP (w/ chemical distributions)
%
%--------------------------------------------------------------------------

%Files pre-loaded: data, plutonic, volcanic, sedimentary
%Outputs: Stats includes median, ranges, number of data points
%age_div = [0:200:4000];

%Plot only mafic volcanics for Nb/La plot?
%agebin = age_hp_box2(data(data.sio2<60 & strcmpi(data.rock_origin,{'volcanic'}),:),[0:200:4000],'All data','sio2')

%data in table format
fig = figure();
%SiO2 value to designate as difference between mafic and felsic
div_value = 60;
textwidth = 16.99757;

for i = 1:length(age_div)-1
    ind = find(age_div(i) <= data.avg_age & data.avg_age < age_div(i+1));
    agebin.n(i) = length(ind);
    agebin.ind{i} = ind;
    agebin.mu(i) = mean(log(data.heat_production(ind)));
    agebin.sigma(i) = std(log(data.heat_production(ind)));
    
    avg_age{i} = data.avg_age(ind);
    heat_production{i} = data.heat_production(ind);
    field_val{i} = data{ind,'sio2'};
end

% heat production versus age
%---------------------------
% linear scale
subplot(8,1,1:2);

[agebin.Qage,agebin.Qhp] = whisker(avg_age,heat_production,'NoPlot');

%Color = [0,0,0]+alpha for greys, [0,0,0] is black, alpha up to 1
sqwavefill(agebin.Qhp,agebin.Qage(:,3),age_div,[0,0,0]);

fprintf('Heat production averages and ages:\n')
[agebin.Qhp(:,3),agebin.Qage(:,3)/1000]

a1 = pettitt(log10(agebin.Qhp(:,3)));
fprintf('Pettitt test (linear-space), age: %f Ma,  p-value: %f\n',age_div(a1(1)+1),a1(3));

a1 = pettitt(log10(agebin.Qhp(:,3)));
fprintf('Pettitt test (log-space), age: %f Ma,  p-value: %f\n',age_div(a1(1)+1),a1(3));
%[h1,p1] = kstest2(data.heat_production(data.avg_age <= age_div(a1(1)+1)), ...
%    data.heat_production(data.avg_age <= age_div(a1(1)+1)));
%[h2,p2] = kstest2(log10(data.heat_production(data.avg_age <= age_div(a1(1)+1))), ...
%    log10(data.heat_production(data.avg_age <= age_div(a1(1)+1))));
%fprintf('K-S test (linear-space), age: %f Ma,  p-value: %f\n',age_div(a1(1)+1),h1);
%fprintf('K-S test (log-space), age: %f Ma,  p-value: %f\n',age_div(a1(1)+1),h2);

%a2 = pettitt(log10(agebin.Qhp(a1(1)+1:end,3)));
%fprintf('Pettitt test (log-space), age: %f Ma,  p-value: %f\n',age_div(a2(1)+1),a2(3));
%[h1,p1] = kstest2(data.heat_production(data.avg_age <= age_div(a2(1)+1)), ...
%    data.heat_production(data.avg_age <= age_div(a2(1)+1)));
%[h2,p2] = kstest2(log10(data.heat_production(data.avg_age <= age_div(a2(1)+1))), ...
%    log10(data.heat_production(data.avg_age <= age_div(a2(1)+1))));
%fprintf('K-S test (linear-space), age: %f Ma,  p-value: %f\n',age_div(a2(1)+1),h1);
%fprintf('K-S test (log-space), age: %f Ma,  p-value: %f\n',age_div(a2(1)+1),h2);

xlim([age_div(1) age_div(end)]);
set(gca,'XTickLabel',{});

set(gca,'Box','on');
ylim([0 ceil(max(agebin.Qhp(:,5)))]);
set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);
%ylabel('Heat Production [\muW m^{-3}]');
ylabel('A [\muW m^{-3}]','FontSize',10)

set(gca,'XAxisLocation','top');
xlabel('Age [Ga]');
set(gca,'XTick',[0:500:4000],'XTickLabel',[0:0.5:4]);

%---------------------------
% log-scale
subplot(8,1,3:4);

%[agebin.Qage,agebin.Qhp] = whisker(avg_age,heat_production,[0.5 0.5 0.5]);
[agebin.Qage,agebin.Qhp] = whisker(avg_age,heat_production,'NoPlot','Scale','log');

%Color = [0,0,0]+alpha for greys, [0,0,0] is black, alpha up to 1
sqwavefill(agebin.Qhp,agebin.Qage(:,3),age_div,[0,0,0])

xlim([age_div(1) age_div(end)]);

set(gca,'Box','on');
set(gca,'XTickLabel',{});
hpax([floor(min(agebin.Qhp(:,1))) ceil(max(agebin.Qhp(:,5)))]);
set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);
%xlabel('Age [Ma]','FontSize',10)
%ylabel('Heat Production [\muW m^{-3}]');
ylabel('A [\muW m^{-3}]','FontSize',10)




% Cloud City plots
%---------------------------
subplot(8,1,5); hold on;
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


% data distribution with age
%---------------------------
% linear scale
subplot(8,1,7);
for i = 1:length(age_div)-1
    N(i) = sum(age_div(i) <= data.avg_age & data.avg_age < age_div(i+1));
end
h = bar(age_div(1:end-1),log10(N),'histc');
hold on;
for i = 1:3
    plot([age_div(1) age_div(end)],[i i],'-', ...
        'LineWidth',0.25,'Color',[0.4 0.4 0.7]);
end
set(h,'FaceColor',[0.7 0.7 0.7])
set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);
hpax([0 4],'y');    
ylabel('N data','FontSize',10);
xlabel('Age [Ma]','FontSize',10);
set(gca,'Box','on')


% Supercontinent cycle
%---------------------------
subplot(8,1,8); hold on;
orogen_hist(0,4000);

% Relative felsic/mafic percentage
%---------------------------
subplot(8,1,6)

for i = 1:length(age_div)-1
    mafic_ind{i} = find(age_div(i) <= data.avg_age & data.avg_age < age_div(i+1) ...
        & data.sio2<=div_value & data.sio2>0);
    felsic_ind{i} = find(age_div(i) <= data.avg_age & data.avg_age < age_div(i+1) ...
        & data.sio2>div_value & data.sio2<=100);
    felsic_num(i) = length(felsic_ind{i});
    mafic_num(i) = length(mafic_ind{i});
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


%pdf
%print('-bestfit',fig,[folder 'hpvsage_raw'],'-dpdf','-r0')
textwidth = 16.99757;
%textwidth = 19;
%textheight = 25;
set(fig, 'units', 'centimeters', 'pos', [0 0 textwidth textwidth])
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters',...
    'PaperSize',[pos(3), pos(4)])

%export_fig '/Users/mgard/Documents/University/PhD/Background Research/Papers/HP vs Age/figures/hpvsage_raw.pdf' -q101

return


% creates cloud city diagrams (mirrored histograms)
function mirrhist(val,N,shift,coord)

if isempty(val)
    warning('Value vector was empty.');
    return;
end

edges = linspace(min(val),max(val),round(sqrt(length(val))));

x = [edges; edges];
x = x(:);

x = [x; flipud(x)];

n = histc(val,edges)';

if N ~=0
    n = n/(N*max(n));
end


y = [0 n(1:end-1); n(1:end-1) 0];
y = y(:);

y = [y; -flipud(y)] + shift;

if strcmp(coord,'y');
    fill(y,x,[0.7 0.7 0.7]);
else
    fill(x,y,[0.7 0.7 0.7]);
end

%This shows lines
% hline = refline([0 65]);
% hline.Color = 'r';
% 
% hline = refline([0 55]);
% hline.Color = 'b';

return