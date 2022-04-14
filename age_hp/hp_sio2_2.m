function [hp_corrected_median] = hp_sio2_2(data,age_div,sio2range,felsic_ind,mafic_ind,folder)

%sio2 correction
cor_val = 75;
hp_corrected = [];
median_hp = [];
no_data = [];

%Median method
x = data.sio2;
y = log10(data.heat_production);

    ind_plut = strcmpi(data.rock_origin,'plutonic');
    ind_volc = strcmpi(data.rock_origin,'volcanic');

r2 = corrcoef(x,y).^2;
fprintf('\nr^2 value for raw HP to SiO2: %f\n',r2(1,2));  
    
for i=1:(length(sio2range)-1)
    ind = (x>=sio2range(i) & x<sio2range(i+1));
    ind_plut_i = ind&ind_plut;
    ind_volc_i = ind&ind_volc;
    median_hp(i) = nanmedian(y(ind));
    no_data(i) = sum(ind);
    no_data_plut(i) = sum(ind_plut_i);
    no_data_volc(i) = sum(ind_volc_i);
end
median_sio2 = (sio2range(2:end)+sio2range(1:end-1))./(sio2range(2:end)-sio2range(1:end-1));
weights_data = no_data ./ (sum(no_data));




x1 = median_sio2(:);
X = [ones(size(x1)) x1];
y = median_hp(:);
weights_data = weights_data(:);

newX = X;
for i = 1:size(X,2)
    newX(:,i) = weights_data.*X(:,i);
end
newY = y.*weights_data;



x = (newX)\(newY);
%First is y int, next is slope





%Make another hp_corrected
hp_corrected_median = [];
for i = 1:size(data.sio2,1)
    hp_corrected_median(i,:) = data.heat_production(i,:).*(10^((cor_val - data.sio2(i,:)).*x(2)));
end


figure()
subplot(2,1,1)
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
    avg_sio2{i} = data.sio2(data.sio2 >= sio2range(i) & data.sio2 < sio2range(i+1) & ~isnan(data.heat_production)& ~isnan(data.sio2),:);
    avg_hp{i} = data.heat_production(data.sio2 >= sio2range(i) & data.sio2 < sio2range(i+1) & ~isnan(data.heat_production)& ~isnan(data.sio2),:);
end
[Qsio2,Qhp] = whisker(avg_sio2,avg_hp,'Color',[0.5 0.5 0.5],'Scale','log');


%Color = [0,0,0]+alpha for greys, [0,0,0] is black, alpha up to 1
plot(0,0)

hold on
%plot(data.sio2,log10(data.heat_production),'.')
sio2range
esio2 = sio2range(1):1:sio2range(end);
ehp = [-2:0.1:2];
n = hist2d(data.sio2,log10(data.heat_production),esio2,ehp);


%imagesc(esio2,ehp,n);
imagesc([min(sio2range)+0.5 max(sio2range)-0.5],[min(ehp)+0.5 max(ehp)-0.5],n)


colormap(flipud(bone));
c = colorbar;
c.Label.String = 'No. Samples';
%bestfitline = plot(data.sio2,x(1)+x(2).*data.sio2,'-r')
hold off
hpax([-2 2]);
xlim([40 85])
golden

sqwavefill(Qhp,Qsio2(:,3),sio2range,[0,0,0])
hold on
plot([44,80],x(1)+x(2).*[44,80],'-r')
hold off
hpax([-2 2]);
xlim([40 85])
golden
ylabel('HP')
xlabel('Age')
set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);
ylabel('A [\muW m^{-3}]','FontSize',10)


x_val = Qsio2(:,3); %x
y_val = Qhp(:,3);   %y
%weights_data %weights
y_val_est = x(1)+x(2).*Qsio2(:,3); %esimated y


%RMSE = sqrt(sum of i = 1 to n of (weight(i).*(x_hat(i)-x(i))^2))

RMSE = 0;
for v = 1:length(Qsio2(:,3))
    RMSE = RMSE + weights_data(v).*((y_val_est(v) - y_val(v)).^2);
end
RMSE = sqrt(RMSE);
fprintf('Root mean square error (log units) = %.2f\n',RMSE)


subplot(2,1,2)
histogram(data.sio2(strcmpi(data.rock_origin,'plutonic')),sio2range,'DisplayStyle','stairs','EdgeColor','b');
hold on
histogram(data.sio2(strcmpi(data.rock_origin,'volcanic')),sio2range,'DisplayStyle','stairs','EdgeColor','r');
histogram(data.sio2,sio2range,'DisplayStyle','stairs','EdgeColor','k');
hold off
ylabel('No. data','FontSize',10);
xlabel('SiO2','FontSize',10);
xlim([40 85])
set(gca,'Box','on')
golden









%lh = legend(bestfitline,['$\mathrm{\log_{10} A = ',num2str(x(2)),'\times SiO_{2}',num2str(x(1)),'}$']);
%set(lh,'Interpreter','latex')
% set(gca,'Units','normalized',...
%         'FontUnits','points',...
%         'FontWeight','normal',...
%         'FontSize',8);
    %lh.FontSize = 9;
    % Get current position (which I used to keep overall size the same)
    %currentLegendPosition = lh.Position;
    % Define new position
    %newLegendPosition = [0.416 0.59 currentLegendPosition([3 4])];
    % Set new position
    %lh.Position = newLegendPosition;


fig = figure();
subplot(4,1,2)
    [agebin.Qage,agebin.Qhp] = whisker(avg_age_felsic,heat_production_felsic,'Color',[0.5 0.5 0.5],'Scale','log');
    [agebin1.Qage,agebin1.Qhp] = whisker(avg_age_mafic,heat_production_mafic,'Color',[0.5 0.5 0.5],'Scale','log');
    hold off
    plot(0,0)
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
%lh2.Position = newLegendPosition;
    
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
    
fig2 = figure;
C = corrcoef(agebin.Qhp(:,3),agebin1.Qhp(:,3));
for i = 1:length(agebin.Qhp(:,3))
    plot([agebin1.Qhp(i,2),agebin1.Qhp(i,4)],[agebin.Qhp(i,3),agebin.Qhp(i,3)],'b-', ...
        [agebin1.Qhp(i,3),agebin1.Qhp(i,3)],[agebin.Qhp(i,2),agebin.Qhp(i,4)],'b-');
    if i == 1
        hold on;
    end
    plot(agebin1.Qhp(i,3),agebin.Qhp(i,3),'ro','MarkerFaceColor','r');
    text(agebin1.Qhp(i,3)+0.05,agebin.Qhp(i,3),num2str(agebin1.Qage(i,3)));
end
m = median([agebin1.Qhp(:,3),agebin.Qhp(:,3)])
plot([-2 2],[-2 2]+m(2) - m(1),'k-');
hpax([-2 1],'x');
hpax([-1 1],'y');
set(gca,'DataAspectRatio',[1 1 1]);
xlabel('Mafic heat production [\muW m^{-3}]');
ylabel('Felsic heat production [\muW m^{-3}]');
text(0,-0.8,['Correlation = ',num2str(C(1,2))]);


figure(fig);
subplot(4,1,4)
    [agebin.Qage,agebin.Qhp] = whisker(avg_age_felsic,heat_production_felsic_corrected,'Color',[0.5 0.5 0.5],'Scale','log');
    [agebin1.Qage,agebin1.Qhp] = whisker(avg_age_mafic,heat_production_mafic_corrected,'Color',[0.5 0.5 0.5],'Scale','log');
    hold off
    plot(0,0)
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


%Heat production distribution of before and after
subplot(4,1,1)
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
%golden;
ylim([0 2000])


felsic_pre_quant = quantile(log10(felsic_pre),[0.16, 0.5, 0.84]);
mafic_pre_quant = quantile(log10(mafic_pre),[0.16, 0.5, 0.84]);

hold on
for i = 1:length(felsic_pre_quant)
    x_val = felsic_pre_quant(i);
    x_range = k(1).BinEdges;
    y_val = k(1).Values;
    for j = 1:(length(x_range)-1)
        if x_val >=  x_range(j) & x_val < x_range(j+1)
            y_peak = y_val(j);
            break;
        else
            continue;
        end
    end
    plot([x_val x_val],[0 y_peak],'-b')
end
hold off

hold on
for i = 1:length(mafic_pre_quant)
    x_val = mafic_pre_quant(i);
    x_range = k(3).BinEdges;
    y_val = k(3).Values;
    for j = 1:(length(x_range)-1)
        if x_val >=  x_range(j) & x_val < x_range(j+1)
            y_peak = y_val(j);
            break;
        else
            continue;
        end
    end
    plot([x_val x_val],[0 y_peak],'-r')
end
hold off

fprintf('Quantiles sio2 (pre):\n')
fprintf('Felsic: ')
for i = 1:length(felsic_pre_quant)
    fprintf('%.3f ',felsic_pre_quant(i)) 
end
fprintf('\n')

fprintf('Mafic: ')
for i = 1:length(mafic_pre_quant)
    fprintf('%.3f ',mafic_pre_quant(i)) 
end
fprintf('\n')


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


subplot(4,1,3)
hold on
k(2) = histogram(log10(felsic_post),'DisplayStyle','stairs','EdgeColor','b');
k(4) = histogram(log10(mafic_post),'DisplayStyle','stairs','EdgeColor','r');
hold off
hpax([-1 2],'x');
ylabel('')
%golden;
ylim([0 2000])


felsic_post_quant = quantile(log10(felsic_post),[0.16, 0.5, 0.84]);
mafic_post_quant = quantile(log10(mafic_post),[0.16, 0.5, 0.84]);

hold on
for i = 1:length(felsic_post_quant)
    x_val = felsic_post_quant(i);
    x_range = k(2).BinEdges;
    y_val = k(2).Values;
    for j = 1:(length(x_range)-1)
        if x_val >=  x_range(j) & x_val < x_range(j+1)
            y_peak = y_val(j);
            break;
        else
            continue;
        end
    end
    plot([x_val x_val],[0 y_peak],'-b')
end
hold off

hold on
for i = 1:length(mafic_post_quant)
    x_val = mafic_post_quant(i);
    x_range = k(4).BinEdges;
    y_val = k(4).Values;
    for j = 1:(length(x_range)-1)
        if x_val >=  x_range(j) & x_val < x_range(j+1)
            y_peak = y_val(j);
            break;
        else
            continue;
        end
    end
    plot([x_val x_val],[0 y_peak],'-r')
end
hold off

fprintf('Quantiles sio2 (post):\n')
fprintf('Felsic: ')
for i = 1:length(felsic_post_quant)
    fprintf('%.3f ',felsic_post_quant(i)) 
end
fprintf('\n')

fprintf('Mafic: ')
for i = 1:length(mafic_post_quant)
    fprintf('%.3f ',mafic_post_quant(i)) 
end
fprintf('\n')
















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
subplot(4,1,2)
    [Qageplut,Qhpplut] = whisker(avg_age_plut,heat_production_plut,'Color',[0.5 0.5 0.5],'Scale','log');
    [Qagevolc,Qhpvolc] = whisker(avg_age_volc,heat_production_volc,'Color',[0.5 0.5 0.5],'Scale','log');
    
    [Qageplut_corrected,Qhpplut_corrected] = whisker(avg_age_plut,heat_production_plut_corrected,[0.5 0.5 0.5],'Color','Scale','log');
    [Qagevolc_corrected,Qhpvolc_corrected] = whisker(avg_age_volc,heat_production_volc_corrected,[0.5 0.5 0.5],'Color','Scale','log');    
    
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


subplot(4,1,4)
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


subplot(4,1,1)


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
%golden;
ylim([0 2000])



plutonic_pre_quant = quantile(log10(plutonic_pre),[0.16, 0.5, 0.84]);
volcanic_pre_quant = quantile(log10(volcanic_pre),[0.16, 0.5, 0.84]);

hold on
for i = 1:length(plutonic_pre_quant)
    x_val = plutonic_pre_quant(i);
    x_range = k(1).BinEdges;
    y_val = k(1).Values;
    for j = 1:(length(x_range)-1)
        if x_val >=  x_range(j) & x_val < x_range(j+1)
            y_peak = y_val(j);
            break;
        else
            continue;
        end
    end
    plot([x_val x_val],[0 y_peak],'-b')
end
hold off

hold on
for i = 1:length(volcanic_pre_quant)
    x_val = volcanic_pre_quant(i);
    x_range = k(3).BinEdges;
    y_val = k(3).Values;
    for j = 1:(length(x_range)-1)
        if x_val >=  x_range(j) & x_val < x_range(j+1)
            y_peak = y_val(j);
            break;
        else
            continue;
        end
    end
    plot([x_val x_val],[0 y_peak],'-r')
end
hold off

fprintf('Quantiles sio2 (pre):\n')
fprintf('Plutonic: ')
for i = 1:length(plutonic_pre_quant)
    fprintf('%.3f ',plutonic_pre_quant(i)) 
end
fprintf('\n')

fprintf('Mafic: ')
for i = 1:length(volcanic_pre_quant)
    fprintf('%.3f ',volcanic_pre_quant(i)) 
end
fprintf('\n')














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







subplot(4,1,3)
hold on
k(2) = histogram(log10(plutonic_post),'DisplayStyle','stairs','EdgeColor','b');
k(4) = histogram(log10(volcanic_post),'DisplayStyle','stairs','EdgeColor','r');
hold off
hpax([-1 2],'x');
ylabel('')
%golden;
ylim([0 2000])

textwidth = 16.99757;
set(gcf, 'units', 'centimeters', 'pos', [0 0 textwidth textwidth])
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters',...
    'PaperSize',[pos(3), pos(4)])




hpax([-1 2],'x');
ylabel('N Data')
%golden;
ylim([0 2000])



plutonic_post_quant = quantile(log10(plutonic_post),[0.16, 0.5, 0.84]);
volcanic_post_quant = quantile(log10(volcanic_post),[0.16, 0.5, 0.84]);

hold on
for i = 1:length(plutonic_post_quant)
    x_val = plutonic_post_quant(i);
    x_range = k(2).BinEdges;
    y_val = k(2).Values;
    for j = 1:(length(x_range)-1)
        if x_val >=  x_range(j) & x_val < x_range(j+1)
            y_peak = y_val(j);
            break;
        else
            continue;
        end
    end
    plot([x_val x_val],[0 y_peak],'-b')
end
hold off

hold on
for i = 1:length(volcanic_post_quant)
    x_val = volcanic_post_quant(i);
    x_range = k(4).BinEdges;
    y_val = k(4).Values;
    for j = 1:(length(x_range)-1)
        if x_val >=  x_range(j) & x_val < x_range(j+1)
            y_peak = y_val(j);
            break;
        else
            continue;
        end
    end
    plot([x_val x_val],[0 y_peak],'-r')
end
hold off

fprintf('Quantiles sio2 (post):\n')
fprintf('Plutonic: ')
for i = 1:length(plutonic_post_quant)
    fprintf('%.3f ',plutonic_post_quant(i)) 
end
fprintf('\n')

fprintf('Mafic: ')
for i = 1:length(volcanic_post_quant)
    fprintf('%.3f ',volcanic_post_quant(i)) 
end
fprintf('\n')







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
    
