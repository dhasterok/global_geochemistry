function new_hp_TA(data,age_div)
    %Selects intervals on the sio2 vs TA plot and does a median calc.
    %These medians are then used to see how the relative abundance of TA
    %rich and poor changes for the time intervals

%Between mafic and felsic rocks - the linear relationship between TA
%and sio2 holds i.e. 45 to 75 sio2
intervals = [floor(min(data.sio2)):1:ceil(max(data.sio2))];

for i = 1:length(intervals)-1
    sio2 = data.sio2(data.sio2 >= intervals(i) & data.sio2 < intervals(i+1),:);
    TA = data.k2o(data.sio2 >= intervals(i) & data.sio2 < intervals(i+1),:)...
        + data.na2o(data.sio2 >= intervals(i) & data.sio2 < intervals(i+1),:);

    median_sio2(i) = median(sio2);
    median_TA(i) = median(TA);
end

figure()
hold on;
% 2-d histogram of tas values
eS = [30:0.5:90];
eTA = [0:0.1:15];
n = hist2d(data.sio2,data.k2o+data.na2o,eS,eTA);
imagesc(eS,eTA,log10(n));
colorbar;
caxis([-0.1 3]);
axis xy;

%Show best fit linear between 45 and 75
plot(median_sio2,median_TA,'-or')
fitvars = polyfit(median_sio2, median_TA, 1);
m = fitvars(1);
c = fitvars(2);
plot(median_sio2,m.*median_sio2 + c,'-g')
hold off
xlim([30 90])
ylim([0 15])


%Now look at time intervals
    %Age div indices
    for i = 1:length(age_div)-1
        ind = find(age_div(i) <= data.avg_age & data.avg_age < age_div(i+1));
        n(i) = length(ind);
        agebin.ind{i} = ind;
        avg_age{i} = data.avg_age(ind);
        heat_production{i} = data.heat_production(ind);
    end
    
    %Loop through ages and find ratio and magnitude ratio of TA to line
    
    number_above = [];
    number_below = [];
    ratio_above = [];
    for i = 1:length(agebin.ind)
        datapoints = data(agebin.ind{i},:);
        avg_age{i} = datapoints.avg_age;
        number_total = (datapoints.na2o + datapoints.k2o)-(m.*datapoints.sio2 + c);
        number_above(i) = length(number_total(number_total>0));
        number_below(i) = length(number_total(number_total<0));
        ratio_above(i) = number_above(i)./(number_above(i) + number_below(i));
        
        datapoints_mafic = datapoints(datapoints.sio2 <= 60,:);
        number_total_mafic = (datapoints_mafic.na2o + datapoints_mafic.k2o)-(m.*datapoints_mafic.sio2 + c);
        number_above_mafic(i) = length(number_total_mafic(number_total_mafic>0));
        number_below_mafic(i) = length(number_total_mafic(number_total_mafic<0));
        
        datapoints_felsic = datapoints(datapoints.sio2 > 60,:);
        number_total_felsic = (datapoints_felsic.na2o + datapoints_felsic.k2o)-(m.*datapoints_felsic.sio2 + c);
        number_above_felsic(i) = length(number_total_felsic(number_total_felsic>0));
        number_below_felsic(i) = length(number_total_felsic(number_total_felsic<0));
        
    end
%     h = bar(age_div(1:end-1),(number_above-number_below)./(number_above+number_below),'histc');
%     set(h,'FaceColor',[0.5 0.5 0.5])
%     ylabel('TA bias [+ve higher TA, -ve lower TA]');
%     xlabel('Age [Ma]');
%     title('TA bias');
%     set(gca,'Box','on')


%     subplot(311)
%     h = bar(age_div(1:end-1),number_above./(number_above + number_below),'histc');
%     set(h,'FaceColor',[0.5 0.5 0.5])
%     ylabel('TA ratio above line');
%     xlabel('Age [Ma]');
%     title('TA bias');
%     set(gca,'Box','on')
%     
%     subplot(312)
%     h = bar(age_div(1:end-1),number_above_felsic./(number_above_felsic + number_below_felsic),'histc');
%     set(h,'FaceColor',[0.5 0.5 0.5])
%     ylabel('TA ratio above line: Felsic');
%     xlabel('Age [Ma]');
%     title('TA bias: Felsic');
%     set(gca,'Box','on')
%     
%     subplot(313)
%     h = bar(age_div(1:end-1),number_above_mafic./(number_above_mafic + number_below_mafic),'histc');
%     set(h,'FaceColor',[0.5 0.5 0.5])
%     ylabel('TA ratio above line: Mafic');
%     xlabel('Age [Ma]');
%     title('TA bias: Mafic');
%     set(gca,'Box','on')


    ratio_TA = (number_above-number_below)./(number_above+number_below);
    ratio_TA_felsic = (number_above_felsic-number_below_felsic)./(number_above_felsic+number_below_felsic);
    ratio_TA_mafic = (number_above_mafic-number_below_mafic)./(number_above_mafic+number_below_mafic);
    
    figure()
    subplot(211)
    plot((age_div(2:end)+age_div(1:end-1))./2,ratio_TA,'-ob')
    title('Total alkali bias: All')
    ylabel('TA bias (+ve above, -ve below)')
    xlabel('Age [Ma]')
    
    subplot(212)
    plot((age_div(2:end)+age_div(1:end-1))./2,ratio_TA_felsic,'-ob')
    hold on
    plot((age_div(2:end)+age_div(1:end-1))./2,ratio_TA_mafic,'-xr')
    hold off
    title('Total alkali bias: Felsic vs Mafic')
    ylabel('TA bias (+ve above, -ve below)')
    xlabel('Age [Ma]')
    
return