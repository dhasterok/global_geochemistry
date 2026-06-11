function liq_age(data,age_div)

for i = 1:length(age_div)-1
    %liquidus indices
    agebin.ind_liq{i} = find(age_div(i) <= data.avg_age & data.avg_age < age_div(i+1) ...
        & ~isnan(data.liquidus));
 
    n_liq(i) = length(agebin.ind_liq{i});

    liq{i} = data.liquidus(agebin.ind_liq{i});
    
    avg_age_liq{i} = data.avg_age(agebin.ind_liq{i});
end

[agebin.Qage_liq,agebin.Qliq] = whisker(avg_age_liq,liq,[0.5 0.5 0.5]);
xlim([age_div(1) age_div(end)]);
ylabel('Liquidus');
title('Liquidus');
set(gca,'Box','on');

str = '';
for i = 1:length(n_liq)
    str = [str,', ',int2str(n_liq(i))];
end
dim = [.125 .6 .3 .3];
annotation('textbox',dim,'String',str,'FitBoxToText','on');


figure()
hold on
g = plot([agebin.Qage_liq(:,3);flipud(agebin.Qage_liq(:,3))]',([agebin.Qliq(:,1);flipud(agebin.Qliq(:,5))])','-');
set(g,'Color',[0.3 0.3 0.3]);
plot(agebin.Qage_liq(:,3),(agebin.Qliq(:,3)),'-ob')  
a(1) = fill([agebin.Qage_liq(:,3);flipud(agebin.Qage_liq(:,3))],([agebin.Qliq(:,2);flipud(agebin.Qliq(:,4))]),'b');
title('Age vs. liquidus')
set(a(1),'facealpha',.3);
set(a(1),'EdgeColor','None');
legend(a(1), 'Variance');
hold off
ylabel('Basalt Liquidus')
xlabel('Age')

return