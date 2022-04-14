function [hppresent, hporigin] = hpinitial3(data,age_div)
%Calculates and plots HP at formation i.e. correction for current day
%values to age it was formed.


Exist_Column = strcmp('avg_age',data.Properties.VariableNames);
val = Exist_Column(Exist_Column==1);
if isempty(val)
    data.avg_age = data.age;
end


%Convert k_ppm to k2o if not listed
ind = (data.k_ppm >= 0 & isnan(data.k2o));
data.k2o(ind,:) = ((2*39.0986 + 15.9994)/(2*39.0986*1e4))*data.k_ppm(ind,:);


%Only plot for k,u, and th positive, avg_age value
ind = find(data.k2o >= 0 & data.u_ppm >= 0 & data.th_ppm >= 0 & data.avg_age >= 0 & data.density_model > 0);


data.hp_present = nan([length(data.sio2),1]);
data.hp_origin = nan([length(data.sio2),1]);
data.hp_ratio = nan([length(data.sio2),1]);

data.hp_ratio = data.hp_corrected./data.heat_production;
u_corrected = data.u_ppm .* data.hp_ratio;
th_corrected = data.th_ppm .* data.hp_ratio;
k_corrected = data.k2o .* data.hp_ratio;


h = waitbar(0,'Please wait...');
for i = 1:length(ind)
    [data.hp_present(ind(i),:),~] = radtime(data.density_model(ind(i),:),k_corrected(ind(i),:),...
        0,0,th_corrected(ind(i),:),u_corrected(ind(i),:),'K2O','Formula','r88');
    [data.hp_origin(ind(i),:),~] = radtime(data.density_model(ind(i),:),k_corrected(ind(i),:),...
        0,0,th_corrected(ind(i),:),u_corrected(ind(i),:),'K2O',...
        'Age',data.avg_age(ind(i),:),'Formula','r88');
    waitbar(i / length(ind))
end

close(h)


%Radtime needs to be vectorised!

% [data.hp_present,~] = radtime(data.density_model,k_corrected,zeros(size(k_corrected)),...
%     zeros(size(k_corrected)),th_corrected,u_corrected,'K2O','Formula','r88');
% [data.hp_origin,~] = radtime(data.density_model,k_corrected,zeros(size(k_corrected)),...
%     zeros(size(k_corrected)),th_corrected,u_corrected,'Age',data.avg_age,'K2O','Age',data.avg_age,'Formula','r88');



for i = 1:length(age_div)-1
       ind = find(age_div(i) <= data.avg_age & data.avg_age < age_div(i+1));
       n(i) = length(ind);
       agebin.ind{i} = ind;
       avg_age{i} = data.avg_age(ind);
       hp_origin{i} = data.hp_origin(ind);
       hp_present{i} = data.hp_present(ind);
end

[agebin.Qage,agebin.Qhp] = whisker(avg_age,hp_origin,[0.5 0.5 0.5],'log');
[agebin1.Qage,agebin1.Qhp] = whisker(avg_age,hp_present,[0.5 0.5 0.5],'log');

out = [agebin.Qage,agebin.Qhp];
save 'decay_adj_age.csv' out -ASCII

subplot(1,2,2)
plot(0,0)
sqwavefill(agebin.Qhp,agebin.Qage(:,3),age_div,[0,0,1],'Original')
xlim([age_div(1) age_div(end)]);
set(gca,'Box','on');
axis square
hpax([-1 2]);

subplot(1,2,1)
plot(0,0)
sqwavefill(agebin1.Qhp,agebin1.Qage(:,3),age_div,[1,0,0],'Present')

xlim([age_div(1) age_div(end)]);

set(gca,'Box','on');
axis square
hpax([-1 2]);

hppresent = data.hp_present;
hporigin = data.hp_origin;

%Corrected in linear space
figure()
title('Corrected in linear space')
plot(0,0)
sqwavefill(10.^(agebin.Qhp),agebin.Qage(:,3),age_div,[0,0,1],'Original')
xlim([age_div(1) age_div(end)]);
ylim([0 10])
set(gca,'Box','on');
golden


return