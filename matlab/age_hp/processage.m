function varargout = processage(data,age_div,varargin)
%   processage(data,age_div) does plots of data at specified age divs
%   example w preloaded data from hp_process (and optional output): 
%   data2 = processage(data,[0:200:4000]);
%
%   data: table format
%   age_div: vector format age divisions

%--------------------------------------------------------------------------
% 1. Preparation and corrections to raw data
%--------------------------------------------------------------------------
% Set background colour to white for figures
set(0,'DefaultFigureColor',[1 1 1])


% Clear the command window and close all plots
close all
clc

% sio2 'felsic' and 'mafic' divisor:
split_fm = 60;
fprintf('Felsic/Mafic split SiO_2: %d wt%%\n',split_fm)

% sio2 correction value - how much sio2 to shfit every sample to
sio2_corr_val = 60;


% Remove Jacqui's data set
ind = strcmpi(data.filename,'ANT_geochem_180824_modified.csv');
data(ind,:) = [];


% Restrict number of samples:
% 1. Heat production value
% 2. An avg_age value
% 3. No oceans
% 4. Igneous samples

ind = (data.heat_production >= 0 & ...
    (~isnan(data.age) | ~isnan(data.age_min) | ~isnan(data.age_max))...
    & ~strcmpi(data.country,'oceanic') & ~strcmpi(data.country,'ocean') & ...
    rockgroup(data,'all igneous'));

%ind2 = strcmpi(data.filename,'Macey_etal2018.csv');

fprintf('Restricting database:\n1. Must have heat production values\n2. Must have an age\n3. No ocean data\n4. Only igneous samples\n')
%fprintf('Shrunk dataset from %d to %d samples\n',size(data,1),size(data(ind,:),1))
%data = data(ind & ~ind2,:);
data = data(ind,:);

% age_var: maximum age uncertainty of points (e.g. 200Ma)
age_var = (age_div(end)-age_div(1))./(length(age_div)-1);
[data.avg_min, data.avg_age, data.age_max] = ...
    age_adjust(data.age_min, data.age, data.age_max, data.age_sd, age_var);

ind = isnan(data.avg_age);
fprintf('No. samples with no avg_age removed: %d\n',length(find(ind)))
data = data(~ind,:);

% Duplicate search
%data = duplicate_check(data);


% Element corrections
data = element_correction(data);

% Silica restriction
sio2range = [40 85];
ind = data.sio2 >= sio2range(1) & data.sio2 <= sio2range(2);
fprintf('\nSilica range restricted from %d to %dwt%%: %d samples removed (%.2f%%)\n',sio2range(1),sio2range(2),size(data,1) - size(data(ind,:),1),(size(data,1) - size(data(ind,:),1))/size(data,1) * 100)
data = data(ind,:);

% Bias check - spatial/file/age
%bias_statistics(data,data.heat_production,age_div)

% Save resulting table as output csv for paper
% Add later - time/memory heavy




%--------------------------------------------------------------------------
% 2. Plots/Functions
%--------------------------------------------------------------------------


% Map plot
%----------
% Sort the data so oldest gets plotted first, youngest last
[~, order] = sort(data.avg_age);
load coastlines

figure()
plot(coastlon,coastlat,'-k')
hold on
scatter(data.longitude(order,:),data.latitude(order,:),10,data.avg_age(order,:),'filled');
hold off
gtspal('contrast.xlsx');
% Map options
ylabel('Latitude','FontSize',10)
xlabel('Longitude','FontSize',10)
xlim([-180 180])
ylim([-90 90])
pbaspect([(1+sqrt(5))/2 1 1])
set(gca,'Box','on');

% Mafic/Felsic correlation: subplot 131
fm_fig = figure;
subplot(131);
fm_corr(data,age_div,split_fm);


% Sample numbers histogram
%--------------------------

fprintf('\n------------------\n')
fprintf('Data distribution\n')
fprintf('------------------\n\n')

fprintf('Age Range     N\n');
fprintf('---------   -----\n');
for i = 1:length(age_div)-1
    N(i) = sum(age_div(i) <= data.avg_age & data.avg_age < age_div(i+1));
    fprintf('%4i %4i  %5i\n',age_div(i),age_div(i+1),N(i));
end
fprintf('-----------------\n');
fprintf('Total:     %6i\n',sum(N));

figure()
h = bar(age_div(1:end-1),log10(N),'histc');
hold on
for i = 1:3
    plot([age_div(1) age_div(end)],[i i],'-', ...
        'LineWidth',0.25,'Color',[0.4 0.4 0.7]);
end
hold off
% Histogram options
set(h,'FaceColor',[0.5 0.5 0.5])
set(gca,'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',8);
hpax([0 4],'y');    
ylabel('N data','FontSize',10);
xlabel('Age [Ma]','FontSize',10);
set(gca,'Box','on')



% Raw data - hp vs age and associated subplots (sio2, etc.)
%----------------------------------------------------------
agebin = hpvsage(data,age_div,split_fm);




% Decay correction plots
%------------------------
data.hp_present = nan([length(data.sio2),1]);
data.hp_origin = nan([length(data.sio2),1]);

data.hpconc_initial = nan([length(data.sio2),1]);
data.hpconc_final = nan([length(data.sio2),1]);

[data.hp_present,data.hp_origin,data.hpconc_initial,data.hpconc_final] = decay_correction(data.avg_age,data.k2o,...
    data.u_ppm,data.th_ppm,data.density_model,age_div);
set(gca,'Units','normalized',...
    'FontUnits','points',...
    'FontWeight','normal',...
    'FontSize',8);
ylabel('A [\muW m^{-3}]','FontSize',10);
xlabel('Age [Ma]','FontSize',10)

subplot(1,2,2)
% Plot median felsic/mafic as dashed lines
hold on
median_felsic = nanmedian(data.hp_origin(data.sio2 > 60 & data.sio2 < 100));
median_mafic = nanmedian(data.hp_origin(data.sio2 <= 60 & data.sio2 > 0));
plot([age_div(1) age_div(end)],log10([median_felsic median_felsic]),'--b')
plot([age_div(1) age_div(end)],log10([median_mafic median_mafic]),'--b')
hold off
fprintf('Median felsic: %f\n',median_felsic)
fprintf('Median mafic: %f\n',median_mafic)




% Mafic/Felsic correlation: subplot 132
tmp = data;
tmp.heat_production = tmp.hp_origin;
figure(fm_fig);
subplot(132);
fm_corr(tmp,age_div,split_fm);
clear tmp


% Rock types plot
%-----------------
% Do top 4 across volcanic/plutonic
fields = rock_hp(data,age_div(1):400:age_div(end));


% SiO2 correction plots
%-----------------------
% Even tighter restriction on sio2 range for sio2 fit
sio2_div = [44:1:80];

% Generate hp_corrected and plot for fit
data.hp_corrected = nan([length(data.sio2),1]);
data.silica_shift = nan([length(data.sio2),1]);
[data.hp_corrected, data.silica_shift] = sio2_correction(data.sio2,...
    data.hp_origin,data.avg_age,age_div,sio2_div,data.rock_origin,sio2_corr_val);


% Generate plots for mafic/felsic and plutonic/volcanic change with sio2
% correction
sio2_correction_plots(data.hp_origin,data.hp_corrected,data.sio2,data.avg_age,...
    data.rock_origin,age_div,split_fm)

% Mafic/Felsic correlation: subplot 133
tmp = data;
tmp.heat_production = tmp.hp_corrected;
figure(fm_fig);
subplot(133);
fm_corr(tmp,age_div,split_fm);
clear tmp

% fft and final trends plot
%---------------------------
sinfit_hp(data.avg_age,data.hp_corrected,data.country,age_div)

% Pettitt test - steps
%a1 = pettitt(log10(agebin.Qhp(:,3)));
%fprintf('Pettitt test (linear-space), age: %f Ma,  p-value: %f\n',age_div(a1(1)+1),a1(3));




% Bias check one last time on adjusted values
bias_statistics(data,data.hp_corrected,age_div)


% Proto aus in vs out
%--------------------

% Prepare binned plots:
avg_age = data.avg_age;
for i = 1:length(age_div)-1
       ind = age_div(i) <= avg_age & avg_age < age_div(i+1);
       ind2 = ind;
       if age_div(i) == 1400 || age_div(i) == 1600 || age_div(i) == 1800
           ind2 = ind & (~strcmpi(data.country,'AU') & ~strcmpi(data.country,'Australia'));
       end
       avg_age_cell{i} = avg_age(ind);
       avg_age_cell_2{i} = avg_age(ind2);
       hpwaus{i} = data.hp_corrected(ind);
       hpnoaus{i} = data.hp_corrected(ind2);
end
figure()
[agebin.Qage,agebin.Qhp] = whisker(avg_age_cell,hpwaus,'color',[0.5 0.5 0.5],'scale','log');
[agebin1.Qage,agebin1.Qhp] = whisker(avg_age_cell_2,hpnoaus,'color',[0.5 0.5 0.5],'scale','log');
subplot(1,2,1)
plot(0,0)
sqwavefill(agebin.Qhp,agebin.Qage(:,3),age_div,[0,0,0],'Original')
title('W aus')
axis square
hpax([-2 2])

subplot(1,2,2)
sqwavefill(agebin1.Qhp,agebin1.Qage(:,3),age_div,[0,0,0],'Original')
title('Without aus (1400-2000)')
axis square
hpax([-2 2])


% Add 'expected' bins
figure()
sqwavefill(agebin1.Qhp,agebin1.Qage(:,3),age_div,[0,0,0],'Original')
x = [100:200:3900];
u_quant = quantile(data.u_ppm(data.avg_age >= 0 & data.avg_age < 200,:),...
    [0.1 0.3 0.50 0.7 0.90]);
th_quant = quantile(data.th_ppm(data.avg_age >= 0 & data.avg_age < 200,:),...
    [0.1 0.3 0.50 0.7 0.90]);
k2o_quant = quantile(data.k2o(data.avg_age >= 0 & data.avg_age < 200,:),...
    [0.1 0.3 0.50 0.7 0.90]);
density_avg = nanmedian(data.density_model(data.avg_age >= 0 & data.avg_age < 200,:));

hp_quant = quantile(data.hp_corrected(data.avg_age >= 0 & data.avg_age < 200,:),[0.05 0.25 0.5 0.75 0.95])
th_u_ratio = nanmedian(data.th_ppm(data.avg_age >= 0 & data.avg_age < 200,:)./data.u_ppm(data.avg_age >= 0 & data.avg_age < 200,:))
k2o_u_ratio = nanmedian(data.k2o(data.avg_age >= 0 & data.avg_age < 200,:)./data.u_ppm(data.avg_age >= 0 & data.avg_age < 200,:))

u_val = (hp_quant .* 10^5)./...
    (density_avg.*((3.48 * 2 * molecularwt('K')/molecularwt('K2O')).*k2o_u_ratio +...
    9.52 + 2.56 .* th_u_ratio))
th_val = u_val .* th_u_ratio
k2o_val = u_val .* k2o_u_ratio


y = [];
for i = 1:length(x)
    [y(i,1),~,~,~] = radtime(density_avg,k2o_val(1),0,0,th_val(1),u_val(1),'Age',x(i),'K2O');
    [y(i,2),~,~,~] = radtime(density_avg,k2o_val(2),0,0,th_val(2),u_val(2),'Age',x(i),'K2O');
    [y(i,3),~,~,~] = radtime(density_avg,k2o_val(3),0,0,th_val(3),u_val(3),'Age',x(i),'K2O');
    [y(i,4),~,~,~] = radtime(density_avg,k2o_val(4),0,0,th_val(4),u_val(4),'Age',x(i),'K2O');
    [y(i,5),~,~,~] = radtime(density_avg,k2o_val(5),0,0,th_val(5),u_val(5),'Age',x(i),'K2O');
end

y

hold on
sqwavefill(log10(y),agebin1.Qage(:,3),age_div,[1,0,0],'Original')
hold off
axis square
hpax([-2 2])

% Strawman model is not true - other processes invariably influence HP
% Either these processes are controlling HP at formation, or the preserved
% rock record is biased to lower values.
% Discuss some influences: How may HP vary in response to some known
% processes
% How may HP be biased to lower values? Thermal stability and selective
% preservation/temperature controls.
% Trying to quantify these models and separate in the observed trends is
% difficult and may encompass future work. Nonetheless, observed HP is
% relatively constant through time, and despite much higher HP
% concentrations in the earth as a whole, rocks formed during this period
% (or preserved in this period) do not exhibit higher HP at formation.

if nargout == 1
    varargout{1} = data;
end

return




% Another figure - expected trend vs observed. Subtract them.
figure
subplot(1,2,1)
plot(0,0)
sqwavefill(agebin1.Qhp,agebin1.Qage(:,3),age_div,[0,0,0],'Original')
title('W aus')
axis square
hpax([-2 2])

% Add trend line - todays (0-200 Ma distribution) and
% project back in steps of 500Ga and fit trend line
initial_HP = nanmedian(data.hp_corrected(data.avg_age >= 0 & data.avg_age < 200,:));
initial_u = nanmedian(data.u_ppm(data.avg_age >= 0 & data.avg_age < 200,:));
initial_th = nanmedian(data.th_ppm(data.avg_age >= 0 & data.avg_age < 200,:));
density_avg = nanmedian(data.density_model(data.avg_age >= 0 & data.avg_age < 200,:));
initial_k = nanmedian(data.k2o(data.avg_age >= 0 & data.avg_age < 200,:));

fprintf('Avg. Density: %f\n',density_avg)
fprintf('Avg. K2O: %f\n',initial_k)
fprintf('Avg. u_ppm: %f\n',initial_u)
fprintf('Avg. th_ppm: %f\n',initial_th)
fprintf('Avg. u235: %f%%\n',0.00711*100)
fprintf('Avg. u238: %f%%\n',0.9928*100)


x = [0:50:4000];
y = [];
for i = 1:length(x)
    [y(i),~,~,~] = radtime(density_avg,initial_k,0,0,initial_th,initial_u,'Age',x(i),'K2O');
end
hold on
plot(x,log10(y),'-k')
%plot(agebin1.Qage(:,3),movmean(agebin1.Qhp(:,3),5),'-k')
hold off
axis square
%hpax([-2 2])


subplot(1,2,2)
size(agebin1.Qhp)
size(y)
size(y(3:4:81))
agebin1.Qhp(:,1)
y(3:4:81)'
agebin1.Qhp(:,1) = agebin1.Qhp(:,1) - log10(y(3:4:81))';
agebin1.Qhp(:,2) = agebin1.Qhp(:,2) - log10(y(3:4:81))';
agebin1.Qhp(:,3) = agebin1.Qhp(:,3) - log10(y(3:4:81))';
agebin1.Qhp(:,4) = agebin1.Qhp(:,4) - log10(y(3:4:81))';
agebin1.Qhp(:,5) = agebin1.Qhp(:,5) - log10(y(3:4:81))';
sqwavefill(agebin1.Qhp,agebin1.Qage(:,3),age_div,[0,0,0],'Original')
axis square



% Plot of map showing spatial bias
geo_age(data,age_div)










% 
% % 
% 
% 
% % Uranium/Thorium through time (shifted for SiO2)
% %-------------------------------------------------
% figure()
% u_ratio = [];
% th_ratio = [];
% for i = 1:length(data.hpconc_initial)
%     temp_th = data.hpconc_initial{i}(:,:,4);
%     temp_u = data.hpconc_initial{i}(:,:,5) + data.hpconc_initial{i}(:,:,6);
%     
%     temp2_th = data.hpconc_final{i}(:,:,4);
%     temp2_u = data.hpconc_final{i}(:,:,5) + data.hpconc_final{i}(:,:,6);
%     
%     u_ratio(i,1) = temp2_u./temp_u;
%     th_ratio(i,1) = temp2_th./temp_th;
% end
% 
% 
% u_corrected = data.u_ppm .* u_ratio;
% th_corrected = data.th_ppm .* th_ratio;
% 
% u_uncorrected = data.u_ppm;
% th_uncorrected = data.th_ppm;
% 
% 
% ratio_uth = th_corrected./u_corrected;
% ratio_uth_uncorrected = th_uncorrected./u_uncorrected;
% 
% % Prepare binned plots:
% avg_age = data.avg_age;
% for i = 1:length(age_div)-1
%        ind = find(age_div(i) <= avg_age & avg_age < age_div(i+1));
%        n(i) = length(ind);
%        agebin.ind{i} = ind;
%        avg_age_cell{i} = avg_age(ind);
%        ratio_uth_cell{i} = ratio_uth(ind);
%        ratio_uth_uncorrected_cell{i} = ratio_uth_uncorrected(ind);
% end
% [agebin.Qage,agebin.Qhp] = whisker(avg_age_cell,ratio_uth_cell,'color',[0.5 0.5 0.5],'scale','log');
% [agebin1.Qage,agebin1.Qhp] = whisker(avg_age_cell,ratio_uth_uncorrected_cell,'color',[0.5 0.5 0.5],'scale','log');
% 
% subplot(1,2,1)
% plot(0,0)
% sqwavefill(agebin1.Qhp,agebin1.Qage(:,3),age_div,[0,0,1],'Original')
% title('Uncorrected')
% 
% subplot(1,2,2)
% sqwavefill(agebin.Qhp,agebin.Qage(:,3),age_div,[1,0,0],'Original')
% title('Corrected')
% 
% 
% if nargout == 1
%     varargout{1} = data;
% end




return