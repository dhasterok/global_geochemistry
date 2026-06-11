function sio2_correction_plots(hp_pre_sio2,hp_post_sio2,sio2,avg_age,rock_origin,age_div,split_fm)
%sio2_correction_plots(hp_pre_sio2,hp_post_sio2,avg_age,split_fm)

fprintf('\n------------------\n')
fprintf('SiO_2 correction plots:\n')
fprintf('------------------\n\n')

warning('off','all')

%--------------------
% Create cell arrays
%--------------------
for i = 1:length(age_div)-1
   avg_age_felsic{i} = avg_age((sio2 > split_fm) & (avg_age >= age_div(i))...
       & (avg_age < age_div(i+1)) & ~isnan(hp_pre_sio2) & ~isnan(hp_post_sio2));
   avg_age_mafic{i} = avg_age((sio2 <= split_fm) & (avg_age >= age_div(i))...
       & (avg_age < age_div(i+1)) & ~isnan(hp_pre_sio2) & ~isnan(hp_post_sio2));
   avg_age_volcanic{i} = avg_age(strcmpi(rock_origin,'volcanic') & (avg_age >= age_div(i))...
       & (avg_age < age_div(i+1)) & ~isnan(hp_pre_sio2) & ~isnan(hp_post_sio2));
   avg_age_plutonic{i} = avg_age(strcmpi(rock_origin,'plutonic') & (avg_age >= age_div(i))...
       & (avg_age < age_div(i+1)) & ~isnan(hp_pre_sio2) & ~isnan(hp_post_sio2));
   
   prehp_felsic{i} = hp_pre_sio2(sio2 > split_fm & avg_age >= age_div(i) ...
       & avg_age < age_div(i+1) & ~isnan(hp_pre_sio2) & ~isnan(hp_post_sio2));
   prehp_mafic{i} = hp_pre_sio2((sio2 <= split_fm) & (avg_age >= age_div(i))...
       & (avg_age < age_div(i+1)) & ~isnan(hp_pre_sio2) & ~isnan(hp_post_sio2));
   prehp_volcanic{i} = hp_pre_sio2(strcmpi(rock_origin,'volcanic') & (avg_age >= age_div(i))...
       & (avg_age < age_div(i+1)) & ~isnan(hp_pre_sio2) & ~isnan(hp_post_sio2));
   prehp_plutonic{i} = hp_pre_sio2(strcmpi(rock_origin,'plutonic') & (avg_age >= age_div(i))...
       & (avg_age < age_div(i+1)) & ~isnan(hp_pre_sio2) & ~isnan(hp_post_sio2));
   
   posthp_felsic{i} = hp_post_sio2(sio2 > split_fm & avg_age >= age_div(i) ...
       & avg_age < age_div(i+1) & ~isnan(hp_pre_sio2) & ~isnan(hp_post_sio2));
   posthp_mafic{i} = hp_post_sio2((sio2 <= split_fm) & (avg_age >= age_div(i))...
       & (avg_age < age_div(i+1)) & ~isnan(hp_pre_sio2) & ~isnan(hp_post_sio2));
   posthp_volcanic{i} = hp_post_sio2(strcmpi(rock_origin,'volcanic') & (avg_age >= age_div(i))...
       & (avg_age < age_div(i+1)) & ~isnan(hp_pre_sio2) & ~isnan(hp_post_sio2));
   posthp_plutonic{i} = hp_post_sio2(strcmpi(rock_origin,'plutonic') & (avg_age >= age_div(i))...
       & (avg_age < age_div(i+1)) & ~isnan(hp_pre_sio2) & ~isnan(hp_post_sio2));
end

mf_plot = figure();
[prefelsic.Qage,prefelsic.Qhp] = whisker(avg_age_felsic,prehp_felsic,'Color',[0.5 0.5 0.5],'Scale','log');
[premafic.Qage,premafic.Qhp] = whisker(avg_age_mafic,prehp_mafic,'Color',[0.5 0.5 0.5],'Scale','log');
[prevolcanic.Qage,prevolcanic.Qhp] = whisker(avg_age_volcanic,prehp_volcanic,'Color',[0.5 0.5 0.5],'Scale','log');
[preplutonic.Qage,preplutonic.Qhp] = whisker(avg_age_plutonic,prehp_plutonic,'Color',[0.5 0.5 0.5],'Scale','log');

[postfelsic.Qage,postfelsic.Qhp] = whisker(avg_age_felsic,posthp_felsic,'Color',[0.5 0.5 0.5],'Scale','log');
[postmafic.Qage,postmafic.Qhp] = whisker(avg_age_mafic,posthp_mafic,'Color',[0.5 0.5 0.5],'Scale','log');
[postvolcanic.Qage,postvolcanic.Qhp] = whisker(avg_age_volcanic,posthp_volcanic,'Color',[0.5 0.5 0.5],'Scale','log');
[postplutonic.Qage,postplutonic.Qhp] = whisker(avg_age_plutonic,posthp_plutonic,'Color',[0.5 0.5 0.5],'Scale','log');


%----------------------
% Hp vs age subplots
%----------------------

% Mafic/Felsic: Pre-sio2 correction age-hp curve
%------------------------------------------------
set(0, 'CurrentFigure', mf_plot)
subplot(4,1,2)
sqwavefill(prefelsic.Qhp,prefelsic.Qage(:,3),age_div,[0,0,1])
hold on
sqwavefill(premafic.Qhp,premafic.Qage(:,3),age_div,[1,0,0])
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
lh2 = legend(h([2 1]),'SiO_{2} > 60','SiO_{2} \leq 60','Interpreter','latex');
set(gca,'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',8);
lh2.FontSize = 9;
currentLegendPosition = lh2.Position;
newLegendPosition = [0.741 0.845 currentLegendPosition([3 4])];
lh2.Position = newLegendPosition;

% Mafic/Felsic: Post-sio2 correction age-hp curve
%-------------------------------------------------
subplot(4,1,4)
sqwavefill(postfelsic.Qhp,postfelsic.Qage(:,3),age_div,[0,0,1])
hold on
sqwavefill(postmafic.Qhp,postmafic.Qage(:,3),age_div,[1,0,0])
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
lh3 = legend(h([2 1]),'SiO_{2} > 60','SiO_{2} \leq 60','Interpreter','latex');
set(gca,'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',8);
lh3.FontSize = 9;
currentLegendPosition = lh3.Position;
newLegendPosition = [0.741 0.246 currentLegendPosition([3 4])];
lh3.Position = newLegendPosition;


% Plutonic/Volcanic: Pre-sio2 correction age-hp curve
%-----------------------------------------------------
pv_plot = figure();
set(0, 'CurrentFigure', pv_plot)
subplot(4,1,2)
sqwavefill(preplutonic.Qhp,preplutonic.Qage(:,3),age_div,[0,0,1])
hold on
sqwavefill(prevolcanic.Qhp,prevolcanic.Qage(:,3),age_div,[1,0,0])
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
lh2 = legend(h([2 1]),'Plutonic','Volcanic','Interpreter','latex');
set(gca,'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',8);
lh2.FontSize = 9;
currentLegendPosition = lh2.Position;
newLegendPosition = [0.741 0.845 currentLegendPosition([3 4])];
lh2.Position = newLegendPosition;

% Plutonic/Volcanic: Post-sio2 correction age-hp curve
%------------------------------------------------------
subplot(4,1,4)
sqwavefill(postplutonic.Qhp,postplutonic.Qage(:,3),age_div,[0,0,1])
hold on
sqwavefill(postvolcanic.Qhp,postvolcanic.Qage(:,3),age_div,[1,0,0])
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
lh3 = legend(h([2 1]),'Plutonic','Volcanic','Interpreter','latex');
set(gca,'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',8);
lh3.FontSize = 9;
currentLegendPosition = lh3.Position;
newLegendPosition = [0.741 0.246 currentLegendPosition([3 4])];
lh3.Position = newLegendPosition;


%------------------------------
% Silica distribution subplots
%------------------------------

hist_felsic_pre = [];
hist_felsic_post = [];
hist_mafic_pre = [];
hist_mafic_post = [];
hist_plutonic_pre = [];
hist_plutonic_post = [];
hist_volcanic_pre = [];
hist_volcanic_post = [];

for i = 1:length(age_div)-1
    hist_felsic_pre = vertcat(hist_felsic_pre,prehp_felsic{i});
    hist_felsic_post = vertcat(hist_felsic_post,posthp_felsic{i});
    hist_mafic_pre = vertcat(hist_mafic_pre,prehp_mafic{i});
    hist_mafic_post = vertcat(hist_mafic_post,posthp_mafic{i});
    
    hist_plutonic_pre = vertcat(hist_plutonic_pre,prehp_plutonic{i});
    hist_plutonic_post = vertcat(hist_plutonic_post,posthp_plutonic{i});
    hist_volcanic_pre = vertcat(hist_volcanic_pre,prehp_volcanic{i});
    hist_volcanic_post = vertcat(hist_volcanic_post,posthp_volcanic{i});
end

% Mafic/Felsic: Pre-sio2 histogram
%----------------------------------
set(0, 'CurrentFigure', mf_plot)
subplot(4,1,1)
hold on
k1 = histogram(log10(hist_felsic_pre),'DisplayStyle','stairs','EdgeColor','b');
k3 = histogram(log10(hist_mafic_pre),'DisplayStyle','stairs','EdgeColor','r');
hold off
hpax([-1 2],'x');
ylabel('N Data')
ylim([0 2000])
% Quantile lines and stats
prefelsic_quant = quantile(log10(hist_felsic_pre),[0.16, 0.5, 0.84]);
premafic_quant = quantile(log10(hist_mafic_pre),[0.16, 0.5, 0.84]);
hold on
y_peak = 0;
for i = 1:length(prefelsic_quant)
    x_val = prefelsic_quant(i);
    x_range = k1.BinEdges;
    y_val = k1.Values;
    for j = 1:(length(x_range)-1)
        if x_val >=  x_range(j) && x_val < x_range(j+1)
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
y_peak = 0;
for i = 1:length(premafic_quant)
    x_val = premafic_quant(i);
    x_range = k3.BinEdges;
    y_val = k3.Values;
    for j = 1:(length(x_range)-1)
        if x_val >=  x_range(j) && x_val < x_range(j+1)
            y_peak = y_val(j);
            break;
        else
            continue;
        end
    end
    plot([x_val x_val],[0 y_peak],'-r')
end
hold off
fprintf('Quantiles (sio2-distribution pre correction):\n')
fprintf('Felsic: ')
for i = 1:length(prefelsic_quant)
    fprintf('%.3f ',prefelsic_quant(i)) 
end
fprintf('\n')

fprintf('Mafic: ')
for i = 1:length(premafic_quant)
    fprintf('%.3f ',premafic_quant(i)) 
end
fprintf('\n')
lh6 = legend('Felsic','Mafic','Interpreter','latex');
set(gca,'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',8);
lh6.FontSize = 9;
currentLegendPosition = lh6.Position;
newLegendPosition = [0.319 0.562 currentLegendPosition([3 4])];
lh6.Position = newLegendPosition;

% Mafic/Felsic: Post-sio2 histogram
%----------------------------------
subplot(4,1,3)
hold on
k2 = histogram(log10(hist_felsic_post),'DisplayStyle','stairs','EdgeColor','b');
k4 = histogram(log10(hist_mafic_post),'DisplayStyle','stairs','EdgeColor','r');
hold off
hpax([-1 2],'x');
ylim([0 2000])

postfelsic_quant = quantile(log10(hist_felsic_post),[0.16, 0.5, 0.84]);
postmafic_quant = quantile(log10(hist_mafic_post),[0.16, 0.5, 0.84]);
hold on
for i = 1:length(postfelsic_quant)
    x_val = postfelsic_quant(i);
    x_range = k2.BinEdges;
    y_val = k2.Values;
    for j = 1:(length(x_range)-1)
        if x_val >=  x_range(j) && x_val < x_range(j+1)
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
for i = 1:length(postmafic_quant)
    x_val = postmafic_quant(i);
    x_range = k4.BinEdges;
    y_val = k4.Values;
    for j = 1:(length(x_range)-1)
        if x_val >=  x_range(j) && x_val < x_range(j+1)
            y_peak = y_val(j);
            break;
        else
            continue;
        end
    end
    plot([x_val x_val],[0 y_peak],'-r')
end
hold off
fprintf('Quantiles (sio2-distribution post correction):\n')
fprintf('Felsic: ')
for i = 1:length(postfelsic_quant)
    fprintf('%.3f ',postfelsic_quant(i)) 
end
fprintf('\n')

fprintf('Mafic: ')
for i = 1:length(postmafic_quant)
    fprintf('%.3f ',postmafic_quant(i)) 
end
fprintf('\n')
lh7 = legend('Felsic','Mafic','Interpreter','latex');
set(gca,'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',8);
lh7.FontSize = 9;
currentLegendPosition = lh7.Position;
newLegendPosition = [0.748 0.574 currentLegendPosition([3 4])];
lh7.Position = newLegendPosition;



% Plutonic/Volcanic: Pre-sio2 histogram
%--------------------------------------
set(0, 'CurrentFigure', pv_plot)
subplot(4,1,1)
hold on
k1 = histogram(log10(hist_plutonic_pre),'DisplayStyle','stairs','EdgeColor','b');
k3 = histogram(log10(hist_volcanic_pre),'DisplayStyle','stairs','EdgeColor','r');
hold off
hpax([-1 2],'x');
ylabel('N Data')
ylim([0 2000])
% Quantile lines and stats
preplutonic_quant = quantile(log10(hist_plutonic_pre),[0.16, 0.5, 0.84]);
prevolcanic_quant = quantile(log10(hist_volcanic_pre),[0.16, 0.5, 0.84]);
hold on
for i = 1:length(preplutonic_quant)
    x_val = preplutonic_quant(i);
    x_range = k1.BinEdges;
    y_val = k1.Values;
    for j = 1:(length(x_range)-1)
        if x_val >=  x_range(j) && x_val < x_range(j+1)
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
for i = 1:length(prevolcanic_quant)
    x_val = prevolcanic_quant(i);
    x_range = k3.BinEdges;
    y_val = k3.Values;
    for j = 1:(length(x_range)-1)
        if x_val >=  x_range(j) && x_val < x_range(j+1)
            y_peak = y_val(j);
            break;
        else
            continue;
        end
    end
    plot([x_val x_val],[0 y_peak],'-r')
end
hold off
fprintf('Quantiles (sio2-distribution pre correction):\n')
fprintf('Plutonic: ')
for i = 1:length(preplutonic_quant)
    fprintf('%.3f ',preplutonic_quant(i)) 
end
fprintf('\n')

fprintf('Volcanic: ')
for i = 1:length(prevolcanic_quant)
    fprintf('%.3f ',prevolcanic_quant(i)) 
end
fprintf('\n')
lh6 = legend('Plutonic','Volcanic','Interpreter','latex');
set(gca,'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',8);
lh6.FontSize = 9;
currentLegendPosition = lh6.Position;
newLegendPosition = [0.319 0.562 currentLegendPosition([3 4])];
lh6.Position = newLegendPosition;

% Plutonic/Volcanic: Post-sio2 histogram
%----------------------------------------
subplot(4,1,3)
hold on
k2 = histogram(log10(hist_plutonic_post),'DisplayStyle','stairs','EdgeColor','b');
k4 = histogram(log10(hist_volcanic_post),'DisplayStyle','stairs','EdgeColor','r');
hold off
hpax([-1 2],'x');
ylim([0 2000])

postplutonic_quant = quantile(log10(hist_plutonic_post),[0.16, 0.5, 0.84]);
postvolcanic_quant = quantile(log10(hist_volcanic_post),[0.16, 0.5, 0.84]);
hold on
for i = 1:length(postplutonic_quant)
    x_val = postplutonic_quant(i);
    x_range = k2.BinEdges;
    y_val = k2.Values;
    for j = 1:(length(x_range)-1)
        if x_val >=  x_range(j) && x_val < x_range(j+1)
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
for i = 1:length(postvolcanic_quant)
    x_val = postvolcanic_quant(i);
    x_range = k4.BinEdges;
    y_val = k4.Values;
    for j = 1:(length(x_range)-1)
        if x_val >=  x_range(j) && x_val < x_range(j+1)
            y_peak = y_val(j);
            break;
        else
            continue;
        end
    end
    plot([x_val x_val],[0 y_peak],'-r')
end
hold off
fprintf('Quantiles (sio2-distribution post correction):\n')
fprintf('Plutonic: ')
for i = 1:length(postplutonic_quant)
    fprintf('%.3f ',postplutonic_quant(i)) 
end
fprintf('\n')

fprintf('Volcanic: ')
for i = 1:length(postvolcanic_quant)
    fprintf('%.3f ',postvolcanic_quant(i)) 
end
fprintf('\n')
lh7 = legend('Plut. - shift','Volc. - shift','Interpreter','latex');
set(gca,'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',8);
lh7.FontSize = 9;
currentLegendPosition = lh7.Position;
newLegendPosition = [0.748 0.574 currentLegendPosition([3 4])];
lh7.Position = newLegendPosition;


%---------------------------------------------
% Combined data - pre vs. post sio2 correction
%---------------------------------------------

% Prepare binned plots:
for i = 1:length(age_div)-1
       ind = find(age_div(i) <= avg_age & avg_age < age_div(i+1));
       n(i) = length(ind);
       agebin.ind{i} = ind;
       avg_age_cell{i} = avg_age(ind);
       hp_pre_sio2_cell{i} = hp_pre_sio2(ind);
       hp_post_sio2_cell{i} = hp_post_sio2(ind);
end

figure()
[agebin.Qage,agebin.Qhp] = whisker(avg_age_cell,hp_pre_sio2_cell,'Color',[0.5 0.5 0.5],'Scale','log');
[agebin1.Qage,agebin1.Qhp] = whisker(avg_age_cell,hp_post_sio2_cell,'Color',[0.5 0.5 0.5],'Scale','log');

% HP at formation
subplot(1,2,1)
plot(0,0)
sqwavefill(agebin.Qhp,agebin.Qage(:,3),age_div,[0,0,1],'Pre-sio2 correction')
xlim([age_div(1) age_div(end)]);
set(gca,'Box','on');
axis square
hpax([-2 2]);

% HP at present day
subplot(1,2,2)
plot(0,0)
sqwavefill(agebin1.Qhp,agebin1.Qage(:,3),age_div,[1,0,0],'Post-sio2 correction')
xlim([age_div(1) age_div(end)]);
set(gca,'Box','on');
axis square
hpax([-2 2]);
title('SiO_2 correction')


warning('on','all')

return