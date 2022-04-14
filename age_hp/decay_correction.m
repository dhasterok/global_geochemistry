function [hp_present,hp_origin,conc_array_initial,conc_array_final] = decay_correction(avg_age,k2o,u_ppm,th_ppm,density_model,age_div)
% decay_correction(age,k2o,u_ppm,th_ppm,density,age_div)
% Calculates and plots HP at given age

%Only plot for k,u, and th positive, avg_age value
ind = find(k2o >= 0 & u_ppm >= 0 & th_ppm >= 0 & avg_age >= 0 & density_model > 0);

% Set up outputs
hp_present = nan([length(k2o),1]);
hp_origin = nan([length(k2o),1]);

u_shift = nan([length(k2o),1]);
th_shift = nan([length(k2o),1]);
conc_array_initial = cell([length(k2o),1]);
conc_array_final = cell([length(k2o),1]);

% Generate present day and age of origin HP
h = waitbar(0,'Please wait...');
for i = 1:length(ind)
    [hp_present(ind(i),:),~,conc_array_initial{ind(i)},~] = radtime(density_model(ind(i),:),k2o(ind(i),:),...
        0,0,th_ppm(ind(i),:),u_ppm(ind(i),:),'K2O','Formula','r88');
    [hp_origin(ind(i),:),~,conc_array_final{ind(i)},~] = radtime(density_model(ind(i),:),k2o(ind(i),:),...
        0,0,th_ppm(ind(i),:),u_ppm(ind(i),:),'K2O',...
        'Age',avg_age(ind(i),:),'Formula','r88');
    if mod(i,100) == 0
        str_decay = sprintf('Decay correction: %d %%',round(100.*(i/length(ind))));
        waitbar(i / length(ind),h,str_decay)
    end
end
delete(h)

% Prepare binned plots:
if nargin < 6
    return
end

for i = 1:length(age_div)-1
       ind = find(age_div(i) <= avg_age & avg_age < age_div(i+1));
       n(i) = length(ind);
       agebin.ind{i} = ind;
       avg_age_cell{i} = avg_age(ind);
       hp_origin_cell{i} = hp_origin(ind);
       hp_present_cell{i} = hp_present(ind);
end

figure()
[agebin.Qage,agebin.Qhp] = whisker(avg_age_cell,hp_origin_cell,'Color',[0.5 0.5 0.5],'Scale','log');
[agebin1.Qage,agebin1.Qhp] = whisker(avg_age_cell,hp_present_cell,'Color',[0.5 0.5 0.5],'Scale','log');

% HP at formation
subplot(1,2,2)
plot(0,0)
sqwavefill(agebin.Qhp,agebin.Qage(:,3),age_div,[0,0,1],'Original')
xlim([age_div(1) age_div(end)]);
set(gca,'Box','on');
axis square
hpax([-2 2]);







% HP at present day
subplot(1,2,1)
plot(0,0)
sqwavefill(agebin1.Qhp,agebin1.Qage(:,3),age_div,[1,0,0],'Present')
xlim([age_div(1) age_div(end)]);
set(gca,'Box','on');
axis square
hpax([-2 2]);
title('Decay correction')


return