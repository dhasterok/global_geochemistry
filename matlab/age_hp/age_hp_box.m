function agebin = age_hp_box(data,age_div,varargin)
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
name = '';

%data in table format
fig = figure();

flag = 0;
if nargin == 3
    flag = 1;
    field = varargin{1};
end

for i = 1:length(age_div)-1
    ind = find(age_div(i) <= data.avg_age & data.avg_age < age_div(i+1));
    n(i) = length(ind);
    agebin.ind{i} = ind;
    avg_age{i} = data.avg_age(ind);
    heat_production{i} = data.heat_production(ind);
    
    if flag
        field_val{i} = data{ind,lower(field)};
    end
end
if flag
    subplot(311);
else
    subplot(211);
end
%[agebin.Qage,agebin.Qhp] = whisker(avg_age,heat_production,'Color',[0.5 0.5 0.5]);
[agebin.Qage,agebin.Qhp] = whisker(avg_age,heat_production,'Color',[0.5 0.5 0.5],'Scale','log');

%Color = [0,0,0]+alpha for greys, [0,0,0] is black, alpha up to 1
sqwavefill(agebin.Qhp,agebin.Qage(:,3),age_div,[0,0,0])

xlim([age_div(1) age_div(end)]);
%ylabel('Heat Production [\muW m^{-3}]');
ylabel('test')
title(name);
set(gca,'Box','on');
hpax([floor(min(agebin.Qhp(:,1))) ceil(max(agebin.Qhp(:,5)))]);
set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);
%xlabel('Age [Ma]','FontSize',10)
ylabel('A [\muW m^{-3}]','FontSize',10)


if flag
    subplot(312);
else
    subplot(212);
end


%Single histogram
    n(n == 0) = NaN;
    logn = log10([n';NaN]);
    % h = bar(age_div,logn,'histc');
    % hpax([0 ceil(max(logn))]);
    % set(h,'FaceColor',[0.5 0.5 0.5]);
    % xlim([age_div(1) age_div(end)]);
    % ylabel('Number of Samples');

%Stacked bar graph
    for i = 1:length(age_div)-1
        ind2 = find(age_div(i) <= data.avg_age & data.avg_age < age_div(i+1) & (strcmpi(data.country,'US') | strcmpi(data.country,'United States')));
        n2(i)=length(ind2);

        ind3 = find(age_div(i) <= data.avg_age & data.avg_age < age_div(i+1) & (strcmpi(data.country,'AU') | strcmpi(data.country,'Australia')));
        n3(i)=length(ind3);

        ind4 = find(age_div(i) <= data.avg_age & data.avg_age < age_div(i+1) & ((~strcmpi(data.country,'AU') | ~strcmpi(data.country,'Australia')) & (~strcmpi(data.country,'US') | ~strcmpi(data.country,'United States'))));
        n4(i)=length(ind4);
    end

    n2(n2 == 0) = NaN;
    logn2 = log10([n2';NaN]);

    n3(n3 == 0) = NaN;
    logn3 = log10([n3';NaN]);

    n4(n4 == 0) = NaN;
    logn4 = log10([n4';NaN]);

    h = bar(age_div+100,[logn2 logn3 logn4],1,'grouped');
    colormap gray

    legend('US','AU','Misc.')
    xlim([age_div(1) age_div(end)]);
    hpax([0 ceil(max(logn))]);
    set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);
    ylabel('No. Samples','FontSize',10);


if ~flag
    xlabel('Age [Ma]','FontSize',10);
    return
end


subplot(313); hold on;
for i = 1:length(field_val)
    
    if isempty(field_val{i})
        continue;
    end
    
    N = 2/(age_div(i+1) - age_div(i));
    %mirrhist(field_val{i},N,agebin.Qage(i,3),'y');
    mirrhist(field_val{i},N,(age_div(i+1) + age_div(i))/2,'y');
end
hold off
xlim([age_div(1) age_div(end)]);
%ylim([30 90])
%ylabel(field);
set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);
xlabel('Age [Ma]','FontSize',10);
ylabel('SiO_2 [wt %]','FontSize',10)
set(gca,'Box','on');

textwidth = 16.99757;
set(fig, 'units', 'centimeters', 'pos', [0 0 textwidth*0.75 textwidth])
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters',...
    'PaperSize',[pos(3), pos(4)])
%pdf
%print(fig,['/Users/mgard/Documents/University/PhD/Background Research/Papers/HP vs Age/figures/' ...
%        'fig' int2str(4)],'-dpdf','-r0')

return


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
    fill(y,x,[0.5 0.5 0.5]);
else
    fill(x,y,[0.5 0.5 0.5]);
end

%This shows lines
% hline = refline([0 65]);
% hline.Color = 'r';
% 
% hline = refline([0 55]);
% hline.Color = 'b';

return