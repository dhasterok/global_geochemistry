function perc_hp_TA(data,age_div)
    %Selects intervals on the sio2 vs TA plot and does a median calc.
    %These medians are then used to see how the relative abundance of TA
    %rich and poor changes for the time intervals

%1. Create the quantile models for each sio2 bin
%------------------------------------------------

%Specify intervals of sio2
intervals = floor(min(data.sio2)):1:ceil(max(data.sio2));
quantile_set = 0.01:0.01:0.99;
quantile_model_x = min(quantile_set):0.001:max(quantile_set);
quantile_model_y = {};

%For each interval of sio2, get the quantiles at 0.01 to 99.99 in steps of
%0.01 and interp between them with spline

plot_quant = [];

for i = 1:length(intervals)-1
    TA = data.k2o(data.sio2 >= intervals(i) & data.sio2 < intervals(i+1),:)...
        + data.na2o(data.sio2 >= intervals(i) & data.sio2 < intervals(i+1),:);
    temp_quant = quantile(TA,quantile_set);
    quantile_model_y{i} = interp1(quantile_set,temp_quant,quantile_model_x,'spline');
    plot_quant(i,:) = quantile(TA,[0.05 0.25 0.5 0.75 0.95]);
end

%Plot quantiles
fig = figure()
TA = data.k2o + data.na2o;
S = data.sio2;
hold on;
% 2-d histogram of tas values
eS = [0:0.5:100];
eTA = [-6:0.2:30];
n = hist2d(S,TA,eS,eTA);
imagesc(eS,eTA,log10(n));
colormap(flipud(bone));
c = colorbar;
c.Label.String = 'No. Samples';

caxis([-0.1 3]);


% plots TAS points from geochemistry
%plot(S,TA,'.')

xlim([39 83]);
ylim([0 14]);
set(gca,'Box','on');
golden;
axis xy;


% % load TAS and UH-Mg fields
% [tasgons,uhmgons] = load_tasgons;
% [ntas,~] = size(tasgons);
% [nuhm,~] = size(uhmgons);
% for i = 1:ntas
%     x = tasgons{i,1}(:,1);
%     y = tasgons{i,1}(:,2);
%     plot(x,y,'w-');
%     %fill(x,y,[0.7 0.7 0.7]);
%     %text(mean(x),mean(y),tasgons{i,3});
% end
% % plots uhmg fields
% for i = 1:nuhm
%     x = uhmgons{i,1}(:,1);
%     y = -uhmgons{i,1}(:,2);
%     plot(x,y,'w-');
%     fill(x,y,[0.7 0.7 0.7]);
%     %text(mean(x),mean(y),uhmgons{i,3});
% end


x = (0.5 + floor(min(data.sio2))):1:(ceil(max(data.sio2)) - 0.5);
plot(x,plot_quant(:,1),'-r','LineWidth',1)
plot(x,plot_quant(:,2),'-r','LineWidth',1)
plot(x,plot_quant(:,3),'-r','LineWidth',1)
plot(x,plot_quant(:,4),'-r','LineWidth',1)
plot(x,plot_quant(:,5),'-r','LineWidth',1)


hold off

set(gca,'Box','on');
set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);
xlabel('SiO_2 [wt%]','FontSize',10)
ylabel('K_2O + Na_2O [wt%]','FontSize',10)
golden
textwidth = 16.99757;
    set(fig, 'units', 'centimeters', 'pos', [0 0 textwidth textwidth/1.6])
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters',...
    'PaperSize',[pos(3), pos(4)])
%export_fig '/Users/mgard/Documents/University/PhD/Background Research/Papers/HP vs Age/figures/tasio2.pdf' -q101



%2. Select the age bins
%------------------------
agebin = {};
cellarrayquant = {};
quantile_age = [];

for i = 1:length(age_div)-1
    agebin{i} = find(age_div(i) <= data.avg_age & data.avg_age < age_div(i+1));
    TA_age = data.k2o(agebin{i},:) + data.na2o(agebin{i},:);
    sio2_age = data.sio2(agebin{i},:);
    [quantile_age(i),cellarrayquant{i}] = med_quant(sio2_age,TA_age,quantile_model_x,quantile_model_y,intervals);
end

fig = figure()
%Make this a square wave plot
%plot((age_div(2:end)+age_div(1:end-1))./2,quantile_age,'-or')

%x = age div = [0 200 200 400 400 600 ... 3600 3800 3800 4000]
%y = quantile age at [100 300 ... 3900]
%age_div
%age = (age_div(2:end)+age_div(1:end-1))./2;
age = (age_div(2:end)+age_div(1:end-1))./2;
lx = length(age);
X = zeros(2*lx,1);
X(1:2:2*lx-1) = age_div(1:lx);
X(2:2:2*lx-2) = age_div(2:lx);
X(end)=age_div(end);

ly = length(quantile_age);
Y = zeros(2*ly,1);
Y(1:2:2*ly-1,:) = quantile_age(1:ly);
Y(2:2:2*ly,:)   = quantile_age(1:ly);
hold on
plot((age_div(2:end)+age_div(1:end-1))./2,quantile_age,'Marker','o','Color',[0,0,0],'linestyle','none') 
plot(X,Y,'linestyle','-','Color',[0,0,0])
hold off
set(gca,'Box','on');
set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);
xlabel('Age [Ma]','FontSize',10)
ylabel('TA Median Percentile','FontSize',10)
ylim([0 0.75])
golden

textwidth = 16.99757;
    set(fig, 'units', 'centimeters', 'pos', [0 0 textwidth textwidth/1.6])
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters',...
    'PaperSize',[pos(3), pos(4)])
%export_fig '/Users/mgard/Documents/University/PhD/Background Research/Papers/HP vs Age/figures/tamedian.pdf' -q101







    %pdf
%    print(fig,['/Users/mgard/Documents/University/PhD/Background Research/Papers/HP vs Age/figures/' ...
%        'fig' int2str(9)],'-dpdf','-r0')


% 
% figure()
% %Why is this messed up? Also might try larger bins?
% [Qage,Qquant] = whisker(agebin,cellarrayquant,[0.5 0.5 0.5]);
% sqwavefill(Qquant,Qage,age_div,'b')
% 
% 
% 
% %MAFIC
% %Specify intervals of sio2
% intervals = floor(min(data.sio2)):1:60;
% quantile_set = 0.01:0.01:0.99;
% quantile_model_x = min(quantile_set):0.001:max(quantile_set);
% quantile_model_y = {};
% temp = data;
% data = data(data.sio2<60,:);
% 
% %For each interval of sio2, get the quantiles at 0.01 to 99.99 in steps of
% %0.01 and interp between them with spline
% 
% for i = 1:length(intervals)-1
%     TA = data.k2o(data.sio2 >= intervals(i) & data.sio2 < intervals(i+1),:)...
%         + data.na2o(data.sio2 >= intervals(i) & data.sio2 < intervals(i+1),:);
%     temp_quant = quantile(TA,quantile_set);
%     quantile_model_y{i} = interp1(quantile_set,temp_quant,quantile_model_x,'spline');
% end
% %2. Select the age bins
% %------------------------
% agebin = {};
% cellarrayquant = {};
% quantile_age = [];
% 
% for i = 1:length(age_div)-1
%     agebin{i} = find(age_div(i) <= data.avg_age & data.avg_age < age_div(i+1));
%     TA_age = data.k2o(agebin{i},:) + data.na2o(agebin{i},:);
%     sio2_age = data.sio2(agebin{i},:);
%     [quantile_age(i),cellarrayquant{i}] = med_quant(sio2_age,TA_age,quantile_model_x,quantile_model_y,intervals);
% end
% 
% q_age_mafic = quantile_age;
% 
% %FELSIC
% %Specify intervals of sio2
% data = temp;
% data = data(data.sio2>=60,:);
% intervals = 60:1:ceil(max(data.sio2));
% quantile_set = 0.01:0.01:0.99;
% quantile_model_x = min(quantile_set):0.001:max(quantile_set);
% quantile_model_y = {};
% 
% %For each interval of sio2, get the quantiles at 0.01 to 99.99 in steps of
% %0.01 and interp between them with spline
% 
% for i = 1:length(intervals)-1
%     TA = data.k2o(data.sio2 >= intervals(i) & data.sio2 < intervals(i+1),:)...
%         + data.na2o(data.sio2 >= intervals(i) & data.sio2 < intervals(i+1),:);
%     temp_quant = quantile(TA,quantile_set);
%     quantile_model_y{i} = interp1(quantile_set,temp_quant,quantile_model_x,'spline');
% end
% %2. Select the age bins
% %------------------------
% agebin = {};
% cellarrayquant = {};
% quantile_age = [];
% 
% for i = 1:length(age_div)-1
%     agebin{i} = find(age_div(i) <= data.avg_age & data.avg_age < age_div(i+1));
%     TA_age = data.k2o(agebin{i},:) + data.na2o(agebin{i},:);
%     sio2_age = data.sio2(agebin{i},:);
%     [quantile_age(i),cellarrayquant{i}] = med_quant(sio2_age,TA_age,quantile_model_x,quantile_model_y,intervals);
% end
% 
% q_age_felsic = quantile_age;
% 
% 
% figure()
% hold on
% plot((age_div(2:end)+age_div(1:end-1))./2,q_age_mafic,'-or')
% plot((age_div(2:end)+age_div(1:end-1))./2,q_age_felsic,'-ob')
% hold off


return


function [quantile_age,cellarrayquant] = med_quant(sio2,TA,quantile_model_x,quantile_model_y,intervals)
    %quantile_model_y accuracy is to 0.001 thus round TA
    TA = round(TA,4);
    quantile_set = zeros(length(TA),1);
    for i = 1:length(intervals)-1
        ind = find(sio2 >= intervals(i) & sio2 < intervals(i+1));
        for j = 1:length(ind)
            [c ind2] = min(abs(quantile_model_y{i}-TA(ind(j))));
            quantile_set(ind(j)) = quantile_model_x(ind2);
        end
    end
    b=find(quantile_set~=0);% get the locations of nonzero. 
    quantile_age = nanmedian(quantile_set(b));%median of the non zero entries
%     quantile_age = nanmedian(quantile_set);
    cellarrayquant = quantile_set;
return