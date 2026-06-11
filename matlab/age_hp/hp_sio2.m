function [hp_corrected_median] = hp_sio2(data,age_div,sio2range,felsic_ind,mafic_ind,folder)

%sio2 correction
cor_val = 75;
hp_corrected = [];

%Median method
x = data.sio2;
y = log10(data.heat_production);
%sio2_range = [45:2:79];

for i=1:(length(sio2range)-1)
    ind = (x>=sio2range(i) & x<sio2range(i+1));
    median_hp(i) = nanmedian(y(ind));
end
median_sio2 = (sio2range(2:end)+sio2range(1:end-1))./(sio2range(2:end)-sio2range(1:end-1));

% figure()
% plot(x,y,'.r')
% hold on
% plot(median_sio2,median_hp,'-og')
% hold off

fit_opts =  fitoptions('Method', 'LinearLeastSquares', 'Robust', 'on');
[fit_result, goodness, output] = fit(median_sio2',median_hp', 'poly1', fit_opts);
%0.03152p1 and -1.839p2 the fit result
%0.0307p1 and -1.7979p2 the robust result
sprintf('Goodness fit SiO2 medians line:\n')
goodness

% hold on
% plot(fit_result,median_sio2',median_hp','-b')
% plot(median_sio2',median_sio2'*0.03147-1.835,'-b')
% hold off

%Make another hp_corrected
hp_corrected_median = [];
for i = 1:size(data.sio2,1)
    hp_corrected_median(i,:) = data.heat_production(i,:).*(10^((cor_val - data.sio2(i,:)).*fit_result.p1));
end

fig = figure();
subplot(3,1,2)

[p,stats] = robustfit(data.sio2,log10(data.heat_production));
hold on
%plot(data.sio2,log10(data.heat_production),'.')
esio2 = sio2range;
ehp = [-2:0.05:2];
n = hist2d(data.sio2,log10(data.heat_production),esio2,ehp);
imagesc(esio2,ehp,n);
colormap(flipud(bone));
c = colorbar;
c.Label.String = 'No. Samples';
bestfitline = plot(data.sio2,p(1)+p(2).*data.sio2,'-r')
hold off
hpax([-2 2]);
xlim([30 90])

%Was here


xlabel('SiO_2 [wt%]','FontSize',10)
ylabel('A [\muW m^{-3}]','FontSize',10)

%And here


for i = 1:length(age_div)-1
        avg_age_felsic{i} = data.avg_age(felsic_ind{i});
        heat_production_felsic{i} = data.heat_production(felsic_ind{i});
        heat_production_felsic_corrected{i} = hp_corrected_median(felsic_ind{i});
        
        avg_age_mafic{i} = data.avg_age(mafic_ind{i});
        heat_production_mafic{i} = data.heat_production(mafic_ind{i});
        heat_production_mafic_corrected{i} = hp_corrected_median(mafic_ind{i});
end


%Calculate whiskers first. Want the next plot to over ride these.
for i = 1:(length(sio2range)-1)
    avg_sio2{i} = data.sio2(data.sio2 >= sio2range(i) & data.sio2 < sio2range(i+1) & ~isnan(data.heat_production),:);
    avg_hp{i} = data.heat_production(data.sio2 >= sio2range(i) & data.sio2 < sio2range(i+1) & ~isnan(data.heat_production),:);
end
[Qsio2,Qhp] = whisker(avg_sio2,avg_hp,'Color',[0.5 0.5 0.5],'Scale','log','NoPlot');


%Color = [0,0,0]+alpha for greys, [0,0,0] is black, alpha up to 1
sqwavefill(Qhp,Qsio2(:,3),sio2range,[0,0,0])
hold on
plot(data.sio2,p(1)+p(2).*data.sio2,'-r')
hold off
hpax([-2 2]);
xlim([40 85])

lh = legend(bestfitline,['$\mathrm{\log_{10} A = ',num2str(p(2)),'\times SiO_{2}',num2str(p(1)),'}$']);
set(lh,'Interpreter','latex')
set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);
    lh.FontSize = 9;
    % Get current position (which I used to keep overall size the same)
    currentLegendPosition = lh.Position;
    % Define new position
    newLegendPosition = [0.416 0.59 currentLegendPosition([3 4])];
    % Set new position
    lh.Position = newLegendPosition;



    subplot(3,1,1)
    [agebin.Qage,agebin.Qhp] = whisker(avg_age_felsic,heat_production_felsic,'Color',[0.5 0.5 0.5],'Scale','log','NoPlot');
    [agebin1.Qage,agebin1.Qhp] = whisker(avg_age_mafic,heat_production_mafic,'Color',[0.5 0.5 0.5],'Scale','log','NoPlot');
    
    sqwavefill(agebin.Qhp,agebin.Qage(:,3),age_div,[0,0,1])
    hold on
    sqwavefill(agebin1.Qhp,agebin1.Qage(:,3),age_div,[1,0,0])
    hold off
    xlim([age_div(1) age_div(end)]);
    
    set(gca,'Box','on');
    hpax([-2 2]);
    set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);
    xlabel('Age [Ma]','FontSize',10)
    ylabel('A [\muW m^{-3}]','FontSize',10)
    
    h = findobj(gca,'Type','patch');
%     [lh2,lhicons] = legend(h([2 1]),'SiO_{2} > 60','SiO_{2} \leq 60');
%     patchinlegend = findobj(lhicons,'type','patch');
%     set(patchinlegend,'facea',0.3)
%     set(lh2,'Interpreter','latex')
%     lh2.FontSize = 9;
%     % Get current position (which I used to keep overall size the same)
%     currentLegendPosition = lh2.Position;
%     % Define new position
%     newLegendPosition = [0.709 0.865 currentLegendPosition([3 4])];
%     % Set new position
%     lh2.Position = newLegendPosition;
    
    lh2 = legend(h([2 1]),'SiO_{2} > 60','SiO_{2} \leq 60','Interpreter','latex');
    set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);
    lh2.FontSize = 9;
    % Get current position (which I used to keep overall size the same)
    currentLegendPosition = lh2.Position;
    % Define new position
    newLegendPosition = [0.741 0.845 currentLegendPosition([3 4])];
    % Set new position
    lh2.Position = newLegendPosition;




subplot(3,1,3)
    [agebin.Qage,agebin.Qhp] = whisker(avg_age_felsic,heat_production_felsic_corrected,'Color',[0.5 0.5 0.5],'Scale','log','NoPlot');
    [agebin1.Qage,agebin1.Qhp] = whisker(avg_age_mafic,heat_production_mafic_corrected,'Color',[0.5 0.5 0.5],'Scale','log','NoPlot');

    sqwavefill(agebin.Qhp,agebin.Qage(:,3),age_div,[0,0,1])
    hold on
    sqwavefill(agebin1.Qhp,agebin1.Qage(:,3),age_div,[1,0,0])
    hold off
    xlim([age_div(1) age_div(end)]);
    
    set(gca,'Box','on');
    hpax([-2 2]);
    set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);
    xlabel('Age [Ma]','FontSize',10)
    ylabel('A [\muW m^{-3}]','FontSize',10)
    
    h = findobj(gca,'Type','patch');
%     [lh3,lhicons] = legend(h([2 1]),'SiO_{2} > 60','SiO_{2} \leq 60');
%     patchinlegend = findobj(lhicons,'type','patch');
%     set(patchinlegend,'facea',0.3)
%     set(lh3,'Interpreter','latex')
%     lh3.FontSize = 9;
%     % Get current position (which I used to keep overall size the same)
%     currentLegendPosition = lh3.Position;
%     % Define new position
%     newLegendPosition = [0.709 0.262 currentLegendPosition([3 4])];
%     % Set new position
%     lh3.Position = newLegendPosition;
    lh3 = legend(h([2 1]),'SiO_{2} > 60','SiO_{2} \leq 60','Interpreter','latex');
    set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);
    lh3.FontSize = 9;
    % Get current position (which I used to keep overall size the same)
    currentLegendPosition = lh3.Position;
    % Define new position
    newLegendPosition = [0.741 0.246 currentLegendPosition([3 4])];
    % Set new position
    lh3.Position = newLegendPosition;




textwidth = 16.99757;
set(fig, 'units', 'centimeters', 'pos', [0 0 textwidth textwidth])
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','centimeters',...
    'PaperSize',[pos(3), pos(4)])
%export_fig '/Users/mgard/Documents/University/PhD/Background Research/Papers/HP vs Age/figures/sio2correct.pdf' -q101


figure()
%Heat production distribution of before and after
subplot(1,2,1)
felsic_pre = [];
felsic_post = [];

mafic_pre = [];
mafic_post = [];

for i = 1:length(age_div)-1
    felsic_pre = vertcat(felsic_pre,heat_production_felsic{i});
    felsic_post = vertcat(felsic_post,heat_production_felsic_corrected{i});
    
    mafic_pre = vertcat(mafic_pre,heat_production_mafic{i});
    mafic_post = vertcat(mafic_post,heat_production_mafic_corrected{i});
end

hold on
k(1) = histogram(log10(felsic_pre),'DisplayStyle','stairs','EdgeColor','b');
k(3) = histogram(log10(mafic_pre),'DisplayStyle','stairs','EdgeColor','r');
hold off

hpax([-1 2],'x');
ylabel('N Data')
golden;
ylim([0 2000])

    lh6 = legend('Felsic','Mafic','Interpreter','latex');
    set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);
    lh6.FontSize = 9;
    % Get current position (which I used to keep overall size the same)
    currentLegendPosition = lh6.Position;
    % Define new position
    newLegendPosition = [0.319 0.562 currentLegendPosition([3 4])];
    % Set new position
    lh6.Position = newLegendPosition;


subplot(1,2,2)
hold on
k(2) = histogram(log10(felsic_post),'DisplayStyle','stairs','EdgeColor','b');
k(4) = histogram(log10(mafic_post),'DisplayStyle','stairs','EdgeColor','r');
hold off
hpax([-1 2],'x');
ylabel('')
golden;
ylim([0 2000])

textwidth = 16.99757;
set(gcf, 'units', 'centimeters', 'pos', [0 0 textwidth textwidth])
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters',...
    'PaperSize',[pos(3), pos(4)])


    lh7 = legend('Plut. - shift','Volc. - shift','Interpreter','latex');
    set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);
    lh7.FontSize = 9;
    % Get current position (which I used to keep overall size the same)
    currentLegendPosition = lh7.Position;
    % Define new position
    newLegendPosition = [0.748 0.574 currentLegendPosition([3 4])];
    % Set new position
    lh7.Position = newLegendPosition;


    
    
    
    
    
    
    
    
    
    
    
    
    

%Plutonic vs Volcanic figure
for i = 1:length(age_div)-1
    ind_plut = find(data.avg_age > age_div(i) & data.avg_age <= age_div(i+1)...
        & (strcmpi(data.rock_origin,'plutonic')|strcmpi(data.rock_origin,'metaplutonic')));
    ind_volc = find(data.avg_age > age_div(i) & data.avg_age <= age_div(i+1)...
        & (strcmpi(data.rock_origin,'volcanic')|strcmpi(data.rock_origin,'metavolcanic')));
    
    avg_age_plut{i} = data.avg_age(ind_plut);
    heat_production_plut{i} = data.heat_production(ind_plut);
    heat_production_plut_corrected{i} = hp_corrected_median(ind_plut);
    
    avg_age_volc{i} = data.avg_age(ind_volc);
    heat_production_volc{i} = data.heat_production(ind_volc);
    heat_production_volc_corrected{i} = hp_corrected_median(ind_volc);
        
end

figure()
subplot(3,2,1:2)
    [Qageplut,Qhpplut] = whisker(avg_age_plut,heat_production_plut,'Color',[0.5 0.5 0.5],'Scale','log','NoPlot');
    [Qagevolc,Qhpvolc] = whisker(avg_age_volc,heat_production_volc,'Color',[0.5 0.5 0.5],'Scale','log','NoPlot');
    
    [Qageplut_corrected,Qhpplut_corrected] = whisker(avg_age_plut,heat_production_plut_corrected,'Color',[0.5 0.5 0.5],'Scale','log','NoPlot');
    [Qagevolc_corrected,Qhpvolc_corrected] = whisker(avg_age_volc,heat_production_volc_corrected,'Color',[0.5 0.5 0.5],'Scale','log','NoPlot');    
    
    hold off
    plot(0,0)
    sqwavefill(Qhpplut,Qageplut(:,3),age_div,[0,0,1])
    hold on
    sqwavefill(Qhpvolc,Qagevolc(:,3),age_div,[1,0,0])
    hold off
    xlim([age_div(1) age_div(end)]);
    
    set(gca,'Box','on');
    hpax([-2 2]);
    set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);
    xlabel('Age [Ma]','FontSize',10)
    ylabel('A [\muW m^{-3}]','FontSize',10)
    
    h = findobj(gca,'Type','patch');
    lh4 = legend(h([2 1]),'Plutonic','Volcanic','Interpreter','latex');
    lh4.FontSize = 9;
    % Get current position (which I used to keep overall size the same)
    currentLegendPosition = lh4.Position;
    % Define new position
    newLegendPosition = [0.762 0.86 currentLegendPosition([3 4])];
    % Set new position
    lh4.Position = newLegendPosition;


subplot(3,2,5:6)
hold off
    plot(0,0)
    sqwavefill(Qhpplut_corrected,Qageplut_corrected(:,3),age_div,[0,0,1])
    hold on
    sqwavefill(Qhpvolc_corrected,Qagevolc_corrected(:,3),age_div,[1,0,0])
    hold off
    xlim([age_div(1) age_div(end)]);
    
    set(gca,'Box','on');
    hpax([-2 2]);
    set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);
    xlabel('Age [Ma]','FontSize',10)
    ylabel('A [\muW m^{-3}]','FontSize',10)
    
    h = findobj(gca,'Type','patch');
    %[lh5,lhicons] = legend(h([2 1]),'Plutonic_corrected','Volcanic_corrected');
    lh5 = legend(h([2 1]),'Plutonic - shift','Volcanic - shift','Interpreter','latex');
    lh5.FontSize = 9;
    % Get current position (which I used to keep overall size the same)
    currentLegendPosition = lh5.Position;
    % Define new position
    newLegendPosition = [0.704 0.26 currentLegendPosition([3 4])];
    % Set new position
    lh5.Position = newLegendPosition;


subplot(3,2,3)


plutonic_pre = [];
plutonic_post = [];

volcanic_pre = [];
volcanic_post = [];

for i = 1:length(age_div)-1
    plutonic_pre = vertcat(plutonic_pre,heat_production_plut{i});
    plutonic_post = vertcat(plutonic_post,heat_production_plut_corrected{i});
    
    volcanic_pre = vertcat(volcanic_pre,heat_production_volc{i});
    volcanic_post = vertcat(volcanic_post,heat_production_volc_corrected{i});
end

hold on
k(1) = histogram(log10(plutonic_pre),'DisplayStyle','stairs','EdgeColor','b');
k(3) = histogram(log10(volcanic_pre),'DisplayStyle','stairs','EdgeColor','r');
hold off

hpax([-1 2],'x');
ylabel('N Data')
golden;
ylim([0 2000])

    lh6 = legend('Plutonic','Volcanic','Interpreter','latex');
    set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);
    lh6.FontSize = 9;
    % Get current position (which I used to keep overall size the same)
    currentLegendPosition = lh6.Position;
    % Define new position
    newLegendPosition = [0.319 0.562 currentLegendPosition([3 4])];
    % Set new position
    lh6.Position = newLegendPosition;







subplot(3,2,4)
hold on
k(2) = histogram(log10(plutonic_post),'DisplayStyle','stairs','EdgeColor','b');
k(4) = histogram(log10(volcanic_post),'DisplayStyle','stairs','EdgeColor','r');
hold off
hpax([-1 2],'x');
ylabel('')
golden;
ylim([0 2000])

textwidth = 16.99757;
set(gcf, 'units', 'centimeters', 'pos', [0 0 textwidth textwidth])
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters',...
    'PaperSize',[pos(3), pos(4)])


    lh7 = legend('Plut. - shift','Volc. - shift','Interpreter','latex');
    set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);
    lh7.FontSize = 9;
    % Get current position (which I used to keep overall size the same)
    currentLegendPosition = lh7.Position;
    % Define new position
    newLegendPosition = [0.748 0.574 currentLegendPosition([3 4])];
    % Set new position
    lh7.Position = newLegendPosition;


return
    
