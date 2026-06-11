function contam(data)
    ind = find(~isnan(data.nb_ppm) & data.nb_ppm>=0 & ...
        ~isnan(data.la_ppm) & data.la_ppm>=0 & ... 
        ~isnan(data.heat_production) & ~isnan(data.u_ppm) & ~isnan(data.th_ppm)...
        & ~isnan(data.k_ppm));
    %plot(sort(data.nb_ppm(ind)./data.la_ppm(ind),'ascend'),'.')
    h = histogram(data.nb_ppm(ind)./data.la_ppm(ind),[0:0.1:3]);
    xlim([0 3])
    title('Nb/La ratio')
    hold on
    %hline = refline([0 1]);
    %hline.Color = 'r';
    line([1,1],ylim,'LineWidth',2,'Color','r')
    hold off
    set(h,'FaceColor',[0.5 0.5 0.5]);
    set(gca,'Box','on');
    %ylim([0 4])
    %xlim([0,length(data.heat_production(ind))])
    
    figure()
    %plot(sort(data.th_ppm(ind)./data.u_ppm(ind),'ascend'),'.')
    h = histogram(data.th_ppm(ind)./data.u_ppm(ind),[0:0.5:20]);
    xlim([0 20])
    title('Th/U ratio')
    hold on
    %hline = refline([0 4]);
    %hline.Color = 'r';
    line([4,4],ylim,'LineWidth',2,'Color','r')
    hold off
    set(h,'FaceColor',[0.5 0.5 0.5]);
    set(gca,'Box','on');
    %ylim([0 16])
    %xlim([0,length(data.heat_production(ind))])
    
    figure()
    subplot(1,2,1)
    %plot(sort(data.k_ppm(ind)./data.u_ppm(ind),'ascend'),'.')
    h = histogram(data.k_ppm(ind)./data.u_ppm(ind),[0:0.1*(10^4):5*(10^4)]);
    xlim([0 5*(10^4)])
    title('Kppm/U ratio')
    %ylim([0,0.4*(10^5)])
    %xlim([0,length(data.heat_production(ind))])
    set(h,'FaceColor',[0.5 0.5 0.5]);
    set(gca,'Box','on');
    
    subplot(1,2,2)
    %k2o in wt%
    %ans * 10000 is in ppm
    %94.1906 is molecular weight of k2o
    %k2 component is 78.1966
    %k2 in ppm = 78.1966./94.1906
    
    h = histogram((data.k2o(ind).*10000.*(78.1966./94.1906))./data.u_ppm(ind),[0:0.1*(10^4):5*(10^4)]);
    %xlim([0 5*(10^4)])
    title('K/U ratio')
    %ylim([0,0.4*(10^5)])
    %xlim([0,length(data.heat_production(ind))])
    set(h,'FaceColor',[0.5 0.5 0.5]);
    set(gca,'Box','on');
    xlim([0 5*(10^4)])
    
    figure()
    subplot(121);
    plot(log10(data.la_ppm(ind)./data.th_ppm(ind)),log10(data.th_ppm(ind)./data.u_ppm(ind)),'.')
    hold on;
    plot(log10([1 10 10 1 1]),log10([2 2 10 10 2]),'-');
    title('Th/La vs Th/U')
    hpax([-1 3],'x')
    xlabel('Th/La')
    hpax([-1 3],'y')
    ylabel('Th/U')
    axis square;
    
    subplot(122);
    title('Th vs U')
    plot(log10(data.u_ppm(ind)),log10(data.th_ppm(ind)),'.');
    hold on;
    plot([-2 2],[-2 2]+log10(2),'-');
    plot([-2 2],[-2 2]+log10(10),'-');
    hpax([-1 2],'x');
    xlabel('U (ppm)');
    hpax([-1 2],'y');
    ylabel('Th (ppm)');
    axis square;
return