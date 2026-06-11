close all
clc
dT = 1;
agebins = (0:dT:4000)';


mean_hp = [];
mean_age = [];
for i = 1:(length(agebins)-1)
    % Get age and hp data for age bin
    age = data2.avg_age(data2.avg_age >= agebins(i) & data2.avg_age < agebins(i+1));
    hp = data2.hp_origin(data2.avg_age >= agebins(i) & data2.avg_age < agebins(i+1));
    
    % Convert to log scale
    hp = log10(hp);
    
    % Take mean
    mean_hp(i) = sum(hp)./length(hp);
    mean_age(i) = sum(age)./length(age);
end

% Remove NaN's
mean_age = mean_age(~isnan(mean_hp))';
mean_hp = mean_hp(~isnan(mean_hp))';

mean_hp = detrend(mean_hp,0);


% mean_age = data2.avg_age;
% mean_hp = data2.hp_origin;
% mean_hp = detrend(mean_hp,0);



% Sort the ages/hp

[mean_age,idx] = sort(mean_age);
mean_hp = mean_hp(idx);



[f gof] = fit(mean_age,mean_hp,'sin4')




figure()
subplot(1,3,1)
plot(mean_age',mean_hp','-b')
hold on
plot(f)
hold off
ylim([-2 2])
axis square


dT = 10;
agebins = (0:dT:4000)';


mean_hp = [];
mean_age = [];
for i = 1:(length(agebins)-1)
    % Get age and hp data for age bin
    age = data2.avg_age(data2.avg_age >= agebins(i) & data2.avg_age < agebins(i+1));
    hp = data2.hp_origin(data2.avg_age >= agebins(i) & data2.avg_age < agebins(i+1));
    
    % Convert to log scale
    hp = log10(hp);
    
    % Take mean
    mean_hp(i) = sum(hp)./length(hp);
    mean_age(i) = sum(age)./length(age);
end

% Remove NaN's
mean_age = mean_age(~isnan(mean_hp))';
mean_hp = mean_hp(~isnan(mean_hp))';

mean_hp = detrend(mean_hp,0);



subplot(1,3,2)
plot(mean_age',mean_hp','-b')
ylim([-2 2])
axis square

subplot(1,3,3)
xlim([0 4000])
h1 = plot(f,'-r');
ylim([-2 2])
axis square
set(h1,'LineWidth',2)



figure()
[S,C,age] = orogen_hist(0,2100);
S = detrend(S,0);
C = detrend(C,0);

%subplot(1,2,1)
plot(age,S)
f = fit(age,S,'sin1')
hold on
plot(f)
hold off

% subplot(1,2,2)
% plot(age,C)
% f = fit(age,C,'sin1');
% hold on
% plot(f)
% hold off



