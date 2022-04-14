function hpVhist(sio2,V,A,rho,SiO2)

nr = length(SiO2) - 1;
for i = 1:nr
    ind = find(SiO2(i) <= sio2 & sio2 < SiO2(i+1));

    subplot(nr,3,3*(nr - i) + 1);
    hold on;
    histogram(log10(A(ind)), ...
        'Normalization','pdf','DisplayStyle','stairs');
  
    if i == 1
        set(gca,'XTick',[-3:3],'XTickLabel',10.^[-3:3]);
        xlabel('Heat Production [\muW m^{-3}]');

        set(gca,'Box','off','TickDir','out');
    elseif i == nr
        set(gca,'XAxisLocation','top');
        set(gca,'XTick',[-3:3],'XTickLabel',10.^[-3:3]);
        set(gca,'Box','off','TickDir','out');
    else
        axis off;
    end
    xlim([-2 2]);
    %pbaspect([3 1 1]);

    subplot(nr,3,3*(nr - i) + 2);
    hold on;
    histogram(rho(ind), ...
        'Normalization','pdf','DisplayStyle','stairs');
    if i == 1
        set(gca,'XTick',[2600:200:3400]);
        xlabel('Density [kg m^{-3}]');
        set(gca,'Box','off','TickDir','out');
    elseif i == nr
        set(gca,'XAxisLocation','top');
        set(gca,'XTick',[2600:200:3400]);
        set(gca,'Box','off','TickDir','out');
    else
        axis off;
    end
    xlim([2500 3500]);
    %pbaspect([3 1 1]);

    subplot(nr,3,3*(nr - i) + 3);
    hold on;
    histogram(V(ind), ...
        'Normalization','pdf','DisplayStyle','stairs');
    if i == 1
        xlabel('Estimated P-Velocity [km s^{-1}]');

        set(gca,'Box','off','TickDir','out');
    elseif i == nr
        set(gca,'XAxisLocation','top');
        xlim([SiO2(i) SiO2(i+1)]);
        set(gca,'Box','off','TickDir','out');
    else
        xlim([SiO2(i) SiO2(i+1)]);
        axis off;
    end
    xlim([5.8 8.4]);
    %pbaspect([3 1 1]);
end

%figure;
%histogram(V);
%xlim([min(SiO2) 8.4]);

return
