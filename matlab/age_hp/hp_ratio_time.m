function hp_ratio_time(data,age_div)

ind = find(isnan(data.k_ppm) & ~isnan(data.k2o));
data.k_ppm(ind) = data(ind,:).k2o * (0.8301 * 100000);


for i = 1:length(age_div)-1
    %th and u indices
    agebin.ind_thu{i} = find(age_div(i) <= data.avg_age & data.avg_age < age_div(i+1) ...
        & ~isnan(data.u_ppm) & ~isnan(data.th_ppm) ...
        & (~strcmp(data.country,'unknown')|~strcmp(data.country,'ocean')));
    %k and u indices
    agebin.ind_ku{i} = find(age_div(i) <= data.avg_age & data.avg_age < age_div(i+1) ...
        & ~isnan(data.u_ppm) & ~isnan(data.k_ppm) ...
        & (~strcmp(data.country,'unknown')|~strcmp(data.country,'ocean')));
    
    n_thu(i) = length(agebin.ind_thu{i});
    n_ku(i) = length(agebin.ind_ku{i});

    thu{i} = data.th_ppm(agebin.ind_thu{i})./data.u_ppm(agebin.ind_thu{i});
    ku{i} = data.k_ppm(agebin.ind_ku{i})./data.u_ppm(agebin.ind_ku{i});
    
    avg_age_thu{i} = data.avg_age(agebin.ind_thu{i});
    avg_age_ku{i} = data.avg_age(agebin.ind_ku{i});
    
end

subplot(1,2,1)
[agebin.Qage_thu,agebin.Qthu] = whisker(avg_age_thu,thu,'Color',[0.5 0.5 0.5],'Scale','log');
xlim([age_div(1) age_div(end)]);
ylabel('Th/U ratio');
title('Th/U ratio');
set(gca,'Box','on');

str = '';
for i = 1:length(n_thu)
    str = [str,', ',int2str(n_thu(i))];
end
dim = [.125 .6 .3 .3];
annotation('textbox',dim,'String',str,'FitBoxToText','on');

subplot(1,2,2)
[agebin.Qage_ku,agebin.Qku] = whisker(avg_age_ku,ku,'Color',[0.5 0.5 0.5],'Scale','log');
xlim([age_div(1) age_div(end)]);
ylabel('K/U ratio');
title('K/U ratio');
set(gca,'Box','on');

str = '';
for i = 1:length(n_ku)
    str = [str,', ',int2str(n_ku(i))];
end
dim = [.575 .6 .3 .3];
annotation('textbox',dim,'String',str,'FitBoxToText','on');



return