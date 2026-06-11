function hp_TA(data,age_div)
    %Best fit sio2 vs TA line:
    %Derived from b = robust(x,y), b(1) yint, b(2) slope
    %Calc this on the fly? Probably better
    sio2 = [0:100];
    TA = 0.1776.*sio2 -5.1874;

    %Age div indices
    for i = 1:length(age_div)-1
        ind = find(~isnan(data.na2o) & (data.k2o + data.na2o) < 20 & (data.k2o + data.na2o) >= 0 & age_div(i) <= data.avg_age & data.avg_age < age_div(i+1) ...
            & ~isnan(data.heat_production) & data.sio2 > 25 & (~strcmp(data.country,'unknown')|~strcmp(data.country,'ocean')));
        n(i) = length(ind);
        agebin.ind{i} = ind;
        avg_age{i} = data.avg_age(ind);
        heat_production{i} = data.heat_production(ind);
    end
    
    %Loop through ages and find ratio and magnitude ratio of TA to line
    
    ratio = [];
    weighted_ratio = [];
    for i = 1:length(agebin.ind)
        datapoints = data(agebin.ind{i},:);
        avg_age{i} = datapoints.avg_age;
        magnitude_variance{i} = (datapoints.na2o + datapoints.k2o) - (0.1776.*datapoints.sio2-5.1874);
        ratio(i) = sum(magnitude_variance{i} > 0)./(sum(magnitude_variance{i} > 0) + sum(magnitude_variance{i} < 0));
        weighted_ratio(i) = sum(magnitude_variance{i})./length(magnitude_variance{i});
        median_ratio(i) = median(magnitude_variance{i});
        uranium{i} = data.u_ppm(agebin.ind{i},:);
        thorium{i} = data.th_ppm(agebin.ind{i},:);
        potassium{i} = data.k2o(agebin.ind{i},:);
        sodium{i} = data.na2o(agebin.ind{i},:);
        silica{i} = data.sio2(agebin.ind{i},:);
        total_alkali{i} = data.na2o(agebin.ind{i},:) + data.k2o(agebin.ind{i},:);
        magnesium{i} = data.mgo(agebin.ind{i},:);
    end

    subplot(1,3,1)
    bar((age_div(2:end)+age_div(1:end-1))./2,ratio)
    
    subplot(1,3,2)
    bar((age_div(2:end)+age_div(1:end-1))./2,weighted_ratio+2)
    
    subplot(1,3,3)
    bar((age_div(2:end)+age_div(1:end-1))./2,median_ratio+2)

    figure()
    [whisker_ratio.Qage,whisker_ratio.Qvar] = whisker(avg_age,magnitude_variance,'Color',[0.5 0.5 0.5]);
    title('Age vs. TA distance from line')
    
    figure()
    hold on
    g = plot([whisker_ratio.Qage(:,3);flipud(whisker_ratio.Qage(:,3))]',([whisker_ratio.Qvar(:,1);flipud(whisker_ratio.Qvar(:,5))])','-');
    set(g,'Color',[0.3 0.3 0.3]);
    plot(whisker_ratio.Qage(:,3),(whisker_ratio.Qvar(:,3)),'-ob')  
    a(1) = fill([whisker_ratio.Qage(:,3);flipud(whisker_ratio.Qage(:,3))],([whisker_ratio.Qvar(:,2);flipud(whisker_ratio.Qvar(:,4))]),'b');
    title('Age vs. TA distance from line')
    
    set(a(1),'facealpha',.3);
    set(a(1),'EdgeColor','None');

    legend(a(1), 'Variance');
    hold off
    
    
    figure()
    [da.Qage,da.Qvar]=whisker(avg_age,uranium,'Color',[0.5 0.5 0.5]);
    title('Uranium')
    
    figure()
    [db.Qage,db.Qvar]=whisker(avg_age,thorium,'Color',[0.5 0.5 0.5]);
    title('Thorium')
    
    figure()
    [dc.Qage,dc.Qvar]=whisker(avg_age,potassium,'Color',[0.5 0.5 0.5]);
    title('Potassium')
    
    figure()
    [ec.Qage,ec.Qvar]=whisker(avg_age,sodium,'Color',[0.5 0.5 0.5]);
    title('Sodium')
    
    figure()
    [fc.Qage,fc.Qvar]=whisker(avg_age,silica,'Color',[0.5 0.5 0.5]);
    title('Silica')
    
    %Silica vs TA sliders through time plot
    figure()
    TA_age(data,[0:200:4000])

    
    %1-2% median bins on the SA vs TA plot and then +ve or -ve
    %Count number above and number below - is that flawed?
    
    
    %sio2 vs Na and K ratio relative to silica correction through time
    %Plot a best fit line like TA vs S and then do distance one? Modify
    %distance one? Or above and below
    
    
    
    %liquidus plot through time, hp corrected through time, corrected hp vs
    %liquidus plot through time AS shaded plots
    
    
    
    
return