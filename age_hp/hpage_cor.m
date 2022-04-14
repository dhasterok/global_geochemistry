function hpage_cor(data,age_div,felsic_ind,mafic_ind)
    fig = figure()
    for i = 1:length(age_div)-1
        avg_age{i} = data.avg_age(vertcat(felsic_ind{i},mafic_ind{i}));
        heat_production{i} = data.heat_production(vertcat(felsic_ind{i},mafic_ind{i}));
        
        avg_age_felsic{i} = data.avg_age(felsic_ind{i});
        heat_production_felsic{i} = data.heat_production(felsic_ind{i});
        
        avg_age_mafic{i} = data.avg_age(mafic_ind{i});
        heat_production_mafic{i} = data.heat_production(mafic_ind{i});
    end
    
    
    subplot(3,1,1)
    [agebin.Qage,agebin.Qhp] = whisker(avg_age,heat_production,'Color',[0.5 0.5 0.5],'Scale','log');
    xlim([age_div(1) age_div(end)]);
    set(gca,'Box','on');
    title('HP vs. Age: All data')
    sqwavefill(agebin.Qhp,agebin.Qage(:,3),age_div,[0,0,0])
    xlim([age_div(1) age_div(end)]);
    ylabel('Heat Production [\muW m^{-3}]');
    set(gca,'Box','on');
    title('All data')
    hpax([floor(min(agebin.Qhp(:,1))) ceil(max(agebin.Qhp(:,5)))]);
    
    subplot(3,1,2)
    [agebin.Qage,agebin.Qhp] = whisker(avg_age_felsic,heat_production_felsic,'Color',[0.5 0.5 0.5],'Scale','log');
    xlim([age_div(1) age_div(end)]);
    set(gca,'Box','on');
    title('HP vs. Age: All data')
    sqwavefill(agebin.Qhp,agebin.Qage(:,3),age_div,[0,0,1])
    xlim([age_div(1) age_div(end)]);
    ylabel('Heat Production [\muW m^{-3}]');
    set(gca,'Box','on');
    title('Felsic (>60 wt% SiO_2)')
    hpax([floor(min(agebin.Qhp(:,1))) ceil(max(agebin.Qhp(:,5)))]);
    
    subplot(3,1,3)
    [agebin.Qage,agebin.Qhp] = whisker(avg_age_mafic,heat_production_mafic,'Color',[0.5 0.5 0.5],'Scale','log');
    xlim([age_div(1) age_div(end)]);
    set(gca,'Box','on');
    title('HP vs. Age: All data')
    sqwavefill(agebin.Qhp,agebin.Qage(:,3),age_div,[1,0,0])
    xlim([age_div(1) age_div(end)]);
    ylabel('Heat Production [\muW m^{-3}]');
    set(gca,'Box','on');
    title('Mafic (<60 wt% SiO_2)')
    hpax([floor(min(agebin.Qhp(:,1))) ceil(max(agebin.Qhp(:,5)))]);
    
return