function hpage_cor2(data,age_div,felsic_ind,mafic_ind)
    fig = figure()
    for i = 1:length(age_div)-1
        avg_age{i} = data.avg_age(vertcat(felsic_ind{i},mafic_ind{i}));
        heat_production{i} = data.hp_corrected(vertcat(felsic_ind{i},mafic_ind{i}));
        
        avg_age_felsic{i} = data.avg_age(felsic_ind{i});
        heat_production_felsic{i} = data.hp_corrected(felsic_ind{i});
        
        avg_age_mafic{i} = data.avg_age(mafic_ind{i});
        heat_production_mafic{i} = data.hp_corrected(mafic_ind{i});
    end
    subplot(2,1,1)
    [agebin.Qage,agebin.Qhp] = whisker(avg_age,heat_production,'Color',[0.5 0.5 0.5],'Scale','log');
    xlim([age_div(1) age_div(end)]);
    set(gca,'Box','on');
    sqwavefill(agebin.Qhp,agebin.Qage(:,3),age_div,[0,0,0])
    xlim([age_div(1) age_div(end)]);
    
    set(gca,'Box','on');
    hpax([floor(min(agebin.Qhp(:,1))) ceil(max(agebin.Qhp(:,5)))]);
    
    set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);
    xlabel('Age [Ma]','FontSize',10);
    ylabel('A [\muW m^{-3}]','FontSize',10)
    
    subplot(2,1,2)
    [agebin.Qage,agebin.Qhp] = whisker(avg_age_felsic,heat_production_felsic,'Color',[0.5 0.5 0.5],'Scale','log');
    [agebin1.Qage,agebin1.Qhp] = whisker(avg_age_mafic,heat_production_mafic,'Color',[0.5 0.5 0.5],'Scale','log');
    hold off
    sqwavefill(agebin.Qhp,agebin.Qage(:,3),age_div,[0,0,1])
    hold on
    sqwavefill(agebin1.Qhp,agebin1.Qage(:,3),age_div,[1,0,0])
    hold off
    xlim([age_div(1) age_div(end)]);
    
    set(gca,'Box','on');
    %hpax([floor(min(vertcat(agebin.Qhp(:,1),agebin1.Qhp(:,1)))) ceil(max(vertcat(agebin.Qhp(:,5),agebin1.Qhp(:,5))))]);
    hpax([-1 2]);
    set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);
    xlabel('Age [Ma]','FontSize',10);
    ylabel('A [\muW m^{-3}]','FontSize',10)
    
    h = findobj(gca,'Type','patch');
    [lh,lhicons] = legend(h([2 1]),'SiO_{2} > 60','SiO_{2} \leq 60');
    patchinlegend = findobj(lhicons,'type','patch');
    set(patchinlegend,'facea',0.3)
    set(lh,'Interpreter','latex')
    lh.FontSize = 10;
    
    textwidth = 16.99757;
    set(fig, 'units', 'centimeters', 'pos', [0 0 textwidth textwidth])
    
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters',...
    'PaperSize',[pos(3), pos(4)])
    %pdf
    %print(fig,['/Users/mgard/Documents/University/PhD/Background Research/Papers/HP vs Age/figures/' ...
    %    'fig' int2str(6)],'-dpdf','-r0')
    
    
    
%     figure;
%     %Cumulative sum - continental growth?
%     plot(tmp_age,flipud(cumsum(flipud(tmp_hp)))./sum(tmp_hp),'-');
%     hold on
%     plot(agebin.Qage(:,3),flipud(cumsum(flipud(10.^(agebin.Qhp(:,3)))))./sum(10.^(agebin.Qhp(:,3))),'-b');
%     plot(agebin1.Qage(:,3),flipud(cumsum(flipud(10.^(agebin1.Qhp(:,3)))))./sum(10.^(agebin1.Qhp(:,3))),'-r');
%     xlim([age_div(1) age_div(end)]);
%     xlabel('Age [Ma]');
%     
%     
%     figure()
%     Plot without log scale (log distribution though makes sense to do it)
%     [agebin.Qage,agebin.Qhp] = whisker(avg_age_felsic,heat_production_felsic,'Color',[0.5 0.5 0.5]);
%     [agebin1.Qage,agebin1.Qhp] = whisker(avg_age_mafic,heat_production_mafic,'Color',[0.5 0.5 0.5]);
%     sqwavefill(agebin.Qhp,agebin.Qage(:,3),age_div,'b')
%     hold on
%     sqwavefill(agebin1.Qhp,agebin1.Qage(:,3),age_div,'r')
%     hold off
return