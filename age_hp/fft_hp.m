function fft_hp(data,varargin)
% fft_hp - moving average of data and then fft/sin1 fit main trend

avg_age = data.avg_age(data.avg_age >= 0 & ~isnan(data.hp_origin),:);
country = data.country(data.avg_age >= 0 & ~isnan(data.hp_origin),:);
sio2 = data.sio2(data.avg_age >= 0 & ~isnan(data.hp_origin),:);
hp_origin = log10(data.hp_origin(data.avg_age >= 0 & ~isnan(data.hp_origin),:));

hp_present = log10(data.hp_present(data.avg_age >= 0 & ~isnan(data.hp_present),:));

% Moving median curve (less sensitive to outliers)
window_step = 1;
window_size = 200;

start_point = window_size/2;
end_point = round(max(avg_age)) - window_size/2;

age = (start_point:window_step:end_point)';
hp = zeros(size(age,1),1);
hp_noaus = zeros(size(age,1),1);

hp_mafic = zeros(size(age,1),1);
hp_felsic = zeros(size(age,1),1);

ind = strcmpi(country,'AU') & avg_age >= 1400 & avg_age < 2000;

for i = 1:length(age)
    %hp(i) = nanmean(hp_origin((avg_age >= (age(i)-(window_size/2))) & (avg_age < (age(i) + (window_size/2)))));
    hp_noaus(i) = nanmedian(hp_origin((avg_age >= (age(i)-(window_size/2))) & (avg_age < (age(i) + (window_size/2))) & ~ind));
    hp(i) = nanmedian(hp_origin((avg_age >= (age(i)-(window_size/2))) & (avg_age < (age(i) + (window_size/2)))));

    hp_mafic(i) = nanmedian(hp_origin((avg_age >= (age(i)-(window_size/2))) & (avg_age < (age(i) + (window_size/2))) & sio2 <= 60 & sio2 >= 38));
    hp_felsic(i) = nanmedian(hp_origin((avg_age >= (age(i)-(window_size/2))) & (avg_age < (age(i) + (window_size/2))) & sio2 > 60 & sio2 <= 80));
end

figure()
r = ksr(avg_age,hp_origin,100);
plot(r.x,r.f,'-r')
hold on
plot(age,hp,'.k')
hold off



return

% fft requires evenly spaced data - interp across NaN fields if lost
hp = interp1(age(~isnan(hp)),hp(~isnan(hp)),age);

hp_noaus = interp1(age(~isnan(hp_noaus)),hp_noaus(~isnan(hp_noaus)),age);

% Detrend it
hp = detrend(hp);
hp = detrend(hp,0);

hp_noaus = detrend(hp_noaus);
hp_noaus = detrend(hp_noaus,0);


% Plot
figure()
subplot(2,2,1)
plot(age,hp)
% hold on
% sinfit = fit(age,hp,'sin1')
% plot(sinfit,'-r')
% hold off
axis square

ylim([-0.6 0.6])
y_lims = ylim;
% hold on
% sc = [173 335; 750 1071; 1500 1800; 2450 2720];
% txt = {'Pangea', 'Rodina', 'Nuna', 'Kenorland'};
% for i = 1:size(sc,1)
%     p = patch('vertices', [sc(i,1), y_lims(1); sc(i,1), y_lims(2); sc(i,2), y_lims(2); sc(i,2) y_lims(1)], ...
%           'faces', [1, 2, 3, 4], ...
%           'FaceColor', 'b', ...
%           'FaceAlpha', 0.2,...
%           'EdgeColor','none');
%     th = text(midpt(sc(i,:)),y_lims(2)-0.125*(y_lims(2)-y_lims(1)),txt{i});
%     th.HorizontalAlignment = 'center';
%     th.VerticalAlignment = 'bottom';
% end
% hold off

ylabel('SiO_2 + decay corrected w Aus')
title('log_1_0(HP) vs. Age')


%fft
Fs = 1000;
L = length(hp) * 20;
L = 2^(nextpow2(L));
Y = fft(hp,L);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
subplot(2,2,2)
plot(f,P1*(L/length(hp))) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (1/Ga)')
ylabel('|P1(f)|')
xlim([0 30])
axis square


% Plot largest frequency
[max_amp, idx] = max(P1*(L/length(hp)));
P3 = angle(Y);
P3 = P3(1:L/2+1);


subplot(2,2,1)
hold on
y = max_amp .* sin(2 * pi * f(idx) * age/1000 - P3(idx));
plot(age,y,'-k')
hold off

th2 = text(100,-0.5,sprintf('log_1_0(HP) = %.3fsin(%.3fage[Ga] - %.3f)', max_amp,2 * pi * f(idx),P3(idx)));
    th.HorizontalAlignment = 'center';
    th.VerticalAlignment = 'bottom';

    
    
subplot(2,2,3)
plot(age,hp_noaus)
axis square
ylim([-0.6 0.6])
ylabel('SiO_2 + decay corrected w/ Aus')
title('log_1_0(HP) vs. Age')

Fs = 1000;
L = length(hp_noaus) * 20;
L = 2^(nextpow2(L));
Y = fft(hp_noaus,L);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;

subplot(2,2,4)
plot(f,P1*(L/length(hp))) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (1/Ga)')
ylabel('|P1(f)|')
xlim([0 30])
axis square


% Plot largest frequency
[max_amp, idx] = max(P1*(L/length(hp_noaus)));
P3 = angle(Y);
P3 = P3(1:L/2+1);


subplot(2,2,3)
hold on
y = max_amp .* sin(2 * pi * f(idx) * age/1000 - P3(idx));
plot(age,y,'-k')
hold off

th2 = text(100,-0.5,sprintf('log_1_0(HP) = %.3fsin(%.3fage[Ga] - %.3f)', max_amp,2 * pi * f(idx),P3(idx)));
    th.HorizontalAlignment = 'center';
    th.VerticalAlignment = 'bottom';


    
% Mafic vs Felsic:

% avg_age vector
% sio2 vector
% hp_origin vector

% Different age vector
hp_mafic = interp1(age(~isnan(hp_mafic)),hp_mafic(~isnan(hp_mafic)),age);
hp_felsic = interp1(age(~isnan(hp_felsic)),hp_felsic(~isnan(hp_felsic)),age);

hp_mafic = detrend(hp_mafic);
hp_mafic = detrend(hp_mafic,0);
hp_felsic = detrend(hp_felsic);
hp_felsic = detrend(hp_felsic,0);




figure()
subplot(2,2,1)
plot(age,hp_mafic,'-r')
axis square
ylim([-0.6 0.6])
ylabel('SiO_2 + decay corrected w Aus - Mafic')
title('log_1_0(HP) vs. Age')

subplot(2,2,3)
plot(age,hp_felsic,'-b')
axis square
ylim([-0.6 0.6])
ylabel('SiO_2 + decay corrected w Aus - Felsic')
title('log_1_0(HP) vs. Age')

subplot(2,2,2)
%fft
Fs = 1000;
L = length(hp_mafic) * 20;
L = 2^(nextpow2(L));
Y = fft(hp_mafic,L);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
subplot(2,2,2)
plot(f,P1*(L/length(hp_mafic))) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (1/Ga)')
ylabel('|P1(f)|')
xlim([0 30])
axis square
% Plot largest frequency
[max_amp, idx] = max(P1*(L/length(hp_mafic)));
P3 = angle(Y);
P3 = P3(1:L/2+1);
subplot(2,2,1)
hold on
y = max_amp .* sin(2 * pi * f(idx) * age/1000 - P3(idx));
plot(age,y,'-k')
hold off
th2 = text(100,-0.5,sprintf('log_1_0(HP) = %.3fsin(%.3fage[Ga] - %.3f)', max_amp,2 * pi * f(idx),P3(idx)));
    th.HorizontalAlignment = 'center';
    th.VerticalAlignment = 'bottom';
    
subplot(2,2,4)
%fft
Fs = 1000;
L = length(hp_felsic) * 20;
L = 2^(nextpow2(L));
Y = fft(hp_felsic,L);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
plot(f,P1*(L/length(hp_felsic))) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (1/Ga)')
ylabel('|P1(f)|')
xlim([0 30])
axis square
% Plot largest frequency
[max_amp, idx] = max(P1*(L/length(hp_felsic)));
P3 = angle(Y);
P3 = P3(1:L/2+1);
subplot(2,2,3)
hold on
y = max_amp .* sin(2 * pi * f(idx) * age/1000 - P3(idx));
plot(age,y,'-k')
hold off
th2 = text(100,-0.5,sprintf('log_1_0(HP) = %.3fsin(%.3fage[Ga] - %.3f)', max_amp,2 * pi * f(idx),P3(idx)));
    th.HorizontalAlignment = 'center';
    th.VerticalAlignment = 'bottom';
    
 

return

%%

































% fft as many data sets and see if any correlations/common frequencies

% Geochem database data:
% 0 SiO2 corrected, decay corrected, with Aus
% 1 SiO2 corrected, decay corrected, without Aus 1400-2000 Ma
% 2 Decay corrected, with Aus (All)
% 3 Decay corrected, without Aus (All)
% 4 Decay corrected, with Aus (Mafic)
% 5 Decay corrected, without Aus (Mafic)
% 6 Decay corrected, with Aus (Felsic)
% 7 Decay corrected, without Aus (Felsic)

% Other data sets:
% 8 Collision (Condie)
% 9 Subduction (Condie)
% 10 eNd data
% 11 eHf data

% Non/partial-cyclical age range additional overlays:
% 12 Super continent age ranges
% 13 Sea level transgression/regression
% 14 Mountain Building
% 15 Zircon frequency/U–Pb detrital zircon crystallization ages
% https://www-sciencedirect-com.proxy.library.adelaide.edu.au/science/article/pii/S1342937X13000506?via%3Dihub

% Varargin:
% 'Overlay': {'SC','Sea','Orogen'}
% Year division, max age, min age: dT, max_age, min_age

max_T = 4000;
min_T = 0;
dT = 1;
agebins = (min_T:dT:max_T)';


figure()
%% 0 - SiO2 corrected, decay corrected, with Aus
mean_hp_0 = [];
mid_age = (agebins(2:end)+agebins(1:end-1))./2;

for i = 1:length(mid_age)
    % Get hp data - sio2 and decay corrected previously
    hp = data.hp_origin(data.avg_age >= agebins(i) & data.avg_age < agebins(i+1));
    % Convert to log scale
    hp = log10(hp);
    
    % Take mean
    mean_hp_0(i) = sum(hp)./length(hp);
end

% Find min and max age can interpolate/use for
findFirstLast = ~isnan(mean_hp_0);
lastIdx = find(findFirstLast(:),1,'last');
firstIdx = find(findFirstLast(:),1,'first');

% Modify mean_hp and mid_age to this
mid_age = mid_age(firstIdx:lastIdx);
mean_hp_0 = mean_hp_0(firstIdx:lastIdx);

% fft requires evenly spaced data - interp across NaN fields lost
mean_hp_0 = interp1(mid_age(~isnan(mean_hp_0)),mean_hp_0(~isnan(mean_hp_0)),mid_age);

% Remove mean/detrend
mean_hp_0 = detrend(mean_hp_0);
mean_hp_0 = detrend(mean_hp_0,0);

% Plot raw data in left subplot
subplot(2,2,1)
plot(mid_age,mean_hp_0,'-r')
ylabel('SiO_2 + decay corrected w Aus')
title('log_1_0(HP) vs. Age')


y_lims = ylim;
hold on
sinfit = fit(mid_age,mean_hp_0,'sin1');
tx = text((max(mid_age) + min(mid_age))/2,y_lims(1)+0.125*(y_lims(2)-y_lims(1)),strcat('Period: ',num2str(1/(sinfit.b1/(2*pi)))));
plot(sinfit,'-k')
hold off
b = gca;
legend(b,'off');



hold on
sc = [173 335; 750 1071; 1500 1800; 2450 2720];
txt = {'Pangea', 'Rodina', 'Nuna', 'Kenorland'};
for i = 1:size(sc,1)
    p = patch('vertices', [sc(i,1), y_lims(1); sc(i,1), y_lims(2); sc(i,2), y_lims(2); sc(i,2) y_lims(1)], ...
          'faces', [1, 2, 3, 4], ...
          'FaceColor', 'b', ...
          'FaceAlpha', 0.2,...
          'EdgeColor','none');
    th = text(midpt(sc(i,:)),y_lims(2)-0.125*(y_lims(2)-y_lims(1)),txt{i});
    th.HorizontalAlignment = 'center';
    th.VerticalAlignment = 'bottom';
end
hold off

ylabel('SiO_2 + decay corrected w Aus')
title('log_1_0(HP) vs. Age')



% Do the fft on the detrended data - zero pad to 20x length and then to
% nearest power of 2 (for speed)
pad_length = length(mean_hp_0) * 20;
pad_length = 2^(nextpow2(pad_length));

fft_0 = fft(mean_hp_0,pad_length);
Pyy = fft_0.*conj(fft_0)/pad_length;
f_0 = (1/dT)/pad_length*(0:floor(pad_length/2))';

% Find 3 highest peaks
[Peak, PeakIdx] = findpeaks(Pyy(1:floor(pad_length/2)+1));
[sorted_peak, sorted_peak_idx] = sort(Peak,'descend');
PeakIdx = PeakIdx(sorted_peak_idx);
sorted_peak = sorted_peak(1:3);
PeakIdx = PeakIdx(1:3);

% Plot the fft frequency curve in right plot
subplot(2,2,2)
plot(f_0,Pyy(1:floor(pad_length/2)+1))
for i = 1:length(PeakIdx)
    text(f_0(PeakIdx(i)), sorted_peak(i), sprintf('%6.2f Ma', 1./f_0(PeakIdx(i))));
end
title('Power spectral density')
xlabel('Frequency (1/Ma)')
xlim([0 0.05])



%% 1 - SiO2 corrected, decay corrected, without Aus 1400-2000 Ma
mean_hp_1 = [];
mid_age = (agebins(2:end)+agebins(1:end-1))./2;
ind_aus = strcmpi(data.country,'AU') & data.avg_age < 2000 & data.avg_age > 1400;
ind_aus = ~ind_aus;


for i = 1:length(mid_age)
    % Get hp data - sio2 and decay corrected previously
    hp = data.hp_origin(data.avg_age >= agebins(i) & data.avg_age < agebins(i+1) & ind_aus);
    % Convert to log scale
    hp = log10(hp);
    
    % Take mean
    mean_hp_1(i) = sum(hp)./length(hp);
end

% Find min and max age can interpolate/use for
findFirstLast = ~isnan(mean_hp_1);
lastIdx = find(findFirstLast(:),1,'last');
firstIdx = find(findFirstLast(:),1,'first');

% Modify mean_hp and mid_age to this
mid_age = mid_age(firstIdx:lastIdx);
mean_hp_1 = mean_hp_1(firstIdx:lastIdx);

% fft requires evenly spaced data - interp across NaN fields lost
mean_hp_1 = interp1(mid_age(~isnan(mean_hp_1)),mean_hp_1(~isnan(mean_hp_1)),mid_age);

% Remove mean/detrend
mean_hp_1 = detrend(mean_hp_1);
mean_hp_1 = detrend(mean_hp_1,0);

% Plot raw data in left subplot
subplot(2,2,3)
hold on
plot(mid_age,mean_hp_1,'-r')
hold off

y_lims = ylim;
hold on
sinfit = fit(mid_age,mean_hp_1,'sin1');
tx = text((max(mid_age) + min(mid_age))/2,y_lims(1)+0.125*(y_lims(2)-y_lims(1)),strcat('Period: ',num2str(1/(sinfit.b1/(2*pi)))));
plot(sinfit,'-k')
hold off
b = gca;
legend(b,'off');



hold on
sc = [173 335; 750 1071; 1500 1800; 2450 2720];
y_lims = ylim;
txt = {'Pangea', 'Rodina', 'Nuna', 'Kenorland'};
for i = 1:size(sc,1)
    p = patch('vertices', [sc(i,1), y_lims(1); sc(i,1), y_lims(2); sc(i,2), y_lims(2); sc(i,2) y_lims(1)], ...
          'faces', [1, 2, 3, 4], ...
          'FaceColor', 'b', ...
          'FaceAlpha', 0.2,...
          'EdgeColor','none');
    th = text(midpt(sc(i,:)),y_lims(2)-0.125*(y_lims(2)-y_lims(1)),txt{i});
    th.HorizontalAlignment = 'center';
    th.VerticalAlignment = 'bottom';
end
hold off
ylabel('SiO_2 + decay corrected w/ Aus')
title('log_1_0(HP) vs. Age')

% Do the fft on the detrended data - zero pad to 20x length and then to
% nearest power of 2 (for speed)
pad_length = length(mean_hp_1) * 20;
pad_length = 2^(nextpow2(pad_length));

fft_1 = fft(mean_hp_1,pad_length);
Pyy = fft_1.*conj(fft_1)/pad_length;
f_1 = (1/dT)/pad_length*(0:floor(pad_length/2))';

% Find 3 highest peaks
[Peak, PeakIdx] = findpeaks(Pyy(1:floor(pad_length/2)+1));
[sorted_peak, sorted_peak_idx] = sort(Peak,'descend');
PeakIdx = PeakIdx(sorted_peak_idx);
sorted_peak = sorted_peak(1:3);
PeakIdx = PeakIdx(1:3);

% Plot the fft frequency curve in right plot
subplot(2,2,4)
plot(f_1,Pyy(1:floor(pad_length/2)+1))
for i = 1:length(PeakIdx)
    text(f_1(PeakIdx(i)), sorted_peak(i), sprintf('%6.2f Ma', 1./f_1(PeakIdx(i))));
end
title('Power spectral density')
xlabel('Frequency (1/Ma)')
xlim([0 0.05])

%% 2 Decay corrected, with Aus
figure()
% Do decay correction here ONLY on samples directly. No SiO2 correction.

%Convert k_ppm to k2o if not listed
ind = (data.k_ppm >= 0 & isnan(data.k2o));
data.k2o(ind,:) = ((2*39.0986 + 15.9994)/(2*39.0986*1e4))*data.k_ppm(ind,:);

%Only plot for k,u, and th positive, avg_age value
ind = find(data.k2o >= 0 & data.u_ppm >= 0 & data.th_ppm >= 0 & data.avg_age >= 0 & data.density_model > 0);

data.hp_origin_nosio2 = nan([length(data.sio2),1]);

% Generate origin hp without sio2 correction
h = waitbar(0,'Please wait...');
for i = 1:length(ind)
    [data.hp_origin_nosio2(ind(i),:),~] = radtime(data.density_model(ind(i),:),data.k2o(ind(i),:),...
        0,0,data.th_ppm(ind(i),:),data.u_ppm(ind(i),:),'K2O',...
        'Age',data.avg_age(ind(i),:),'Formula','r88');
    waitbar(i / length(ind))
end
close(h)

% fft
mean_hp_2 = [];
mid_age = (agebins(2:end)+agebins(1:end-1))./2;


for i = 1:length(mid_age)
    % Get hp data - sio2 and decay corrected previously
    hp = data.hp_origin_nosio2(data.avg_age >= agebins(i) & data.avg_age < agebins(i+1));
    % Convert to log scale
    hp = log10(hp);
    
    % Take mean
    mean_hp_2(i) = sum(hp)./length(hp);
end

% Find min and max age can interpolate/use for
findFirstLast = ~isnan(mean_hp_2);
lastIdx = find(findFirstLast(:),1,'last');
firstIdx = find(findFirstLast(:),1,'first');

% Modify mean_hp and mid_age to this
mid_age = mid_age(firstIdx:lastIdx);
mean_hp_2 = mean_hp_2(firstIdx:lastIdx);

% fft requires evenly spaced data - interp across NaN fields lost
mean_hp_2 = interp1(mid_age(~isnan(mean_hp_2)),mean_hp_2(~isnan(mean_hp_2)),mid_age);

% Remove mean/detrend
mean_hp_2 = detrend(mean_hp_2);
mean_hp_2 = detrend(mean_hp_2,0);

% Plot raw data in left subplot
subplot(2,2,1)
hold on
plot(mid_age,mean_hp_2,'-r')
hold off

y_lims = ylim;
hold on
sinfit = fit(mid_age,mean_hp_2,'sin1');
tx = text((max(mid_age) + min(mid_age))/2,y_lims(1)+0.125*(y_lims(2)-y_lims(1)),strcat('Period: ',num2str(1/(sinfit.b1/(2*pi)))));
plot(sinfit,'-k')
hold off
b = gca;
legend(b,'off');



hold on
sc = [173 335; 750 1071; 1500 1800; 2450 2720];
y_lims = ylim;
txt = {'Pangea', 'Rodina', 'Nuna', 'Kenorland'};
for i = 1:size(sc,1)
    p = patch('vertices', [sc(i,1), y_lims(1); sc(i,1), y_lims(2); sc(i,2), y_lims(2); sc(i,2) y_lims(1)], ...
          'faces', [1, 2, 3, 4], ...
          'FaceColor', 'b', ...
          'FaceAlpha', 0.2,...
          'EdgeColor','none');
    th = text(midpt(sc(i,:)),y_lims(2)-0.125*(y_lims(2)-y_lims(1)),txt{i});
    th.HorizontalAlignment = 'center';
    th.VerticalAlignment = 'bottom';
end
hold off
ylabel('SiO_2 + decay corrected w/ Aus')
title('log_1_0(HP) vs. Age')

% Do the fft on the detrended data - zero pad to 20x length and then to
% nearest power of 2 (for speed)
pad_length = length(mean_hp_2) * 20;
pad_length = 2^(nextpow2(pad_length));

fft_2 = fft(mean_hp_2,pad_length);
Pyy = fft_2.*conj(fft_2)/pad_length;
f_2 = (1/dT)/pad_length*(0:floor(pad_length/2))';

% Find 3 highest peaks
[Peak, PeakIdx] = findpeaks(Pyy(1:floor(pad_length/2)+1));
[sorted_peak, sorted_peak_idx] = sort(Peak,'descend');
PeakIdx = PeakIdx(sorted_peak_idx);
sorted_peak = sorted_peak(1:3);
PeakIdx = PeakIdx(1:3);

% Plot the fft frequency curve in right plot
subplot(2,2,2)
plot(f_2,Pyy(1:floor(pad_length/2)+1))
for i = 1:length(PeakIdx)
    text(f_2(PeakIdx(i)), sorted_peak(i), sprintf('%6.2f Ma', 1./f_2(PeakIdx(i))));
end
title('Power spectral density')
xlabel('Frequency (1/Ma)')
xlim([0 0.05])

%% 3 - Decay corrected, without Aus 1400-2000 Ma
mean_hp_3 = [];
mid_age = (agebins(2:end)+agebins(1:end-1))./2;
ind_aus = strcmpi(data.country,'AU') & data.avg_age < 2000 & data.avg_age > 1400;
ind_aus = ~ind_aus;


for i = 1:length(mid_age)
    % Get hp data - sio2 and decay corrected previously
    hp = data.hp_origin_nosio2(data.avg_age >= agebins(i) & data.avg_age < agebins(i+1) & ind_aus);
    % Convert to log scale
    hp = log10(hp);
    
    % Take mean
    mean_hp_3(i) = sum(hp)./length(hp);
end

% Find min and max age can interpolate/use for
findFirstLast = ~isnan(mean_hp_3);
lastIdx = find(findFirstLast(:),1,'last');
firstIdx = find(findFirstLast(:),1,'first');

% Modify mean_hp and mid_age to this
mid_age = mid_age(firstIdx:lastIdx);
mean_hp_3 = mean_hp_3(firstIdx:lastIdx);

% fft requires evenly spaced data - interp across NaN fields lost
mean_hp_3 = interp1(mid_age(~isnan(mean_hp_3)),mean_hp_3(~isnan(mean_hp_3)),mid_age);

% Remove mean/detrend
mean_hp_3 = detrend(mean_hp_3);
mean_hp_3 = detrend(mean_hp_3,0);

% Plot raw data in left subplot
subplot(2,2,3)
hold on
plot(mid_age,mean_hp_3,'-r')
hold off

y_lims = ylim;
hold on
sinfit = fit(mid_age,mean_hp_3,'sin1');
tx = text((max(mid_age) + min(mid_age))/2,y_lims(1)+0.125*(y_lims(2)-y_lims(1)),strcat('Period: ',num2str(1/(sinfit.b1/(2*pi)))));
plot(sinfit,'-k')
hold off
b = gca;
legend(b,'off');


hold on
sc = [173 335; 750 1071; 1500 1800; 2450 2720];
y_lims = ylim;
txt = {'Pangea', 'Rodina', 'Nuna', 'Kenorland'};
for i = 1:size(sc,1)
    p = patch('vertices', [sc(i,1), y_lims(1); sc(i,1), y_lims(2); sc(i,2), y_lims(2); sc(i,2) y_lims(1)], ...
          'faces', [1, 2, 3, 4], ...
          'FaceColor', 'b', ...
          'FaceAlpha', 0.2,...
          'EdgeColor','none');
    th = text(midpt(sc(i,:)),y_lims(2)-0.125*(y_lims(2)-y_lims(1)),txt{i});
    th.HorizontalAlignment = 'center';
    th.VerticalAlignment = 'bottom';
end
hold off
ylabel('Decay corrected w/ Aus')
title('log_1_0(HP) vs. Age')

% Do the fft on the detrended data - zero pad to 20x length and then to
% nearest power of 2 (for speed)
pad_length = length(mean_hp_3) * 20;
pad_length = 2^(nextpow2(pad_length));

fft_3 = fft(mean_hp_3,pad_length);
Pyy = fft_3.*conj(fft_3)/pad_length;
f_3 = (1/dT)/pad_length*(0:floor(pad_length/2))';

% Find 3 highest peaks
[Peak, PeakIdx] = findpeaks(Pyy(1:floor(pad_length/2)+1));
[sorted_peak, sorted_peak_idx] = sort(Peak,'descend');
PeakIdx = PeakIdx(sorted_peak_idx);
sorted_peak = sorted_peak(1:3);
PeakIdx = PeakIdx(1:3);

% Plot the fft frequency curve in right plot
subplot(2,2,4)
plot(f_3,Pyy(1:floor(pad_length/2)+1))
for i = 1:length(PeakIdx)
    text(f_3(PeakIdx(i)), sorted_peak(i), sprintf('%6.2f Ma', 1./f_3(PeakIdx(i))));
end
title('Power spectral density')
xlabel('Frequency (1/Ma)')
xlim([0 0.05])





%% 4 Decay corrected, with Aus (Mafic)
figure()
% fft
mean_hp_4 = [];
mid_age = (agebins(2:end)+agebins(1:end-1))./2;


for i = 1:length(mid_age)
    % Get hp data - sio2 and decay corrected previously
    hp = data.hp_origin_nosio2(data.avg_age >= agebins(i) & data.avg_age < agebins(i+1) & data.sio2 <= 60);
    % Convert to log scale
    hp = log10(hp);
    
    % Take mean
    mean_hp_4(i) = sum(hp)./length(hp);
end

% Find min and max age can interpolate/use for
findFirstLast = ~isnan(mean_hp_4);
lastIdx = find(findFirstLast(:),1,'last');
firstIdx = find(findFirstLast(:),1,'first');

% Modify mean_hp and mid_age to this
mid_age = mid_age(firstIdx:lastIdx);
mean_hp_4 = mean_hp_4(firstIdx:lastIdx);

% fft requires evenly spaced data - interp across NaN fields lost
mean_hp_4 = interp1(mid_age(~isnan(mean_hp_4)),mean_hp_4(~isnan(mean_hp_4)),mid_age);

% Remove mean/detrend
mean_hp_4 = detrend(mean_hp_4);
mean_hp_4 = detrend(mean_hp_4,0);

% Plot raw data in left subplot
subplot(2,2,1)
hold on
plot(mid_age,mean_hp_4,'-r')
hold off

y_lims = ylim;
hold on
sinfit = fit(mid_age,mean_hp_4,'sin1');
tx = text((max(mid_age) + min(mid_age))/2,y_lims(1)+0.125*(y_lims(2)-y_lims(1)),strcat('Period: ',num2str(1/(sinfit.b1/(2*pi)))));
plot(sinfit,'-k')
hold off
b = gca;
legend(b,'off');


hold on
sc = [173 335; 750 1071; 1500 1800; 2450 2720];
y_lims = ylim;
txt = {'Pangea', 'Rodina', 'Nuna', 'Kenorland'};
for i = 1:size(sc,1)
    p = patch('vertices', [sc(i,1), y_lims(1); sc(i,1), y_lims(2); sc(i,2), y_lims(2); sc(i,2) y_lims(1)], ...
          'faces', [1, 2, 3, 4], ...
          'FaceColor', 'b', ...
          'FaceAlpha', 0.2,...
          'EdgeColor','none');
    th = text(midpt(sc(i,:)),y_lims(2)-0.125*(y_lims(2)-y_lims(1)),txt{i});
    th.HorizontalAlignment = 'center';
    th.VerticalAlignment = 'bottom';
end
hold off
ylabel('SiO_2 + decay corrected w/ Aus (Mafic)')
title('log_1_0(HP) vs. Age')

% Do the fft on the detrended data - zero pad to 20x length and then to
% nearest power of 2 (for speed)
pad_length = length(mean_hp_4) * 20;
pad_length = 2^(nextpow2(pad_length));

fft_4 = fft(mean_hp_4,pad_length);
Pyy = fft_4.*conj(fft_4)/pad_length;
f_4 = (1/dT)/pad_length*(0:floor(pad_length/2))';

% Find 3 highest peaks
[Peak, PeakIdx] = findpeaks(Pyy(1:floor(pad_length/2)+1));
[sorted_peak, sorted_peak_idx] = sort(Peak,'descend');
PeakIdx = PeakIdx(sorted_peak_idx);
sorted_peak = sorted_peak(1:3);
PeakIdx = PeakIdx(1:3);

% Plot the fft frequency curve in right plot
subplot(2,2,2)
plot(f_4,Pyy(1:floor(pad_length/2)+1))
for i = 1:length(PeakIdx)
    text(f_4(PeakIdx(i)), sorted_peak(i), sprintf('%6.2f Ma', 1./f_4(PeakIdx(i))));
end
title('Power spectral density')
xlabel('Frequency (1/Ma)')
xlim([0 0.05])

%% 5 - Decay corrected, without Aus 1400-2000 Ma (Mafic)
mean_hp_5 = [];
mid_age = (agebins(2:end)+agebins(1:end-1))./2;
ind_aus = strcmpi(data.country,'AU') & data.avg_age < 2000 & data.avg_age > 1400;
ind_aus = ~ind_aus;


for i = 1:length(mid_age)
    % Get hp data - sio2 and decay corrected previously
    hp = data.hp_origin_nosio2(data.avg_age >= agebins(i) & data.avg_age < agebins(i+1) & ind_aus & data.sio2 <= 60);
    % Convert to log scale
    hp = log10(hp);
    
    % Take mean
    mean_hp_5(i) = sum(hp)./length(hp);
end

% Find min and max age can interpolate/use for
findFirstLast = ~isnan(mean_hp_5);
lastIdx = find(findFirstLast(:),1,'last');
firstIdx = find(findFirstLast(:),1,'first');

% Modify mean_hp and mid_age to this
mid_age = mid_age(firstIdx:lastIdx);
mean_hp_5 = mean_hp_5(firstIdx:lastIdx);

% fft requires evenly spaced data - interp across NaN fields lost
mean_hp_5 = interp1(mid_age(~isnan(mean_hp_5)),mean_hp_5(~isnan(mean_hp_5)),mid_age);

% Remove mean/detrend
mean_hp_5 = detrend(mean_hp_5);
mean_hp_5 = detrend(mean_hp_5,0);

% Plot raw data in left subplot
subplot(2,2,3)
hold on
plot(mid_age,mean_hp_5,'-r')
hold off

y_lims = ylim;
hold on
sinfit = fit(mid_age,mean_hp_5,'sin1');
tx = text((max(mid_age) + min(mid_age))/2,y_lims(1)+0.125*(y_lims(2)-y_lims(1)),strcat('Period: ',num2str(1/(sinfit.b1/(2*pi)))));
plot(sinfit,'-k')
hold off
b = gca;
legend(b,'off');


hold on
sc = [173 335; 750 1071; 1500 1800; 2450 2720];
y_lims = ylim;
txt = {'Pangea', 'Rodina', 'Nuna', 'Kenorland'};
for i = 1:size(sc,1)
    p = patch('vertices', [sc(i,1), y_lims(1); sc(i,1), y_lims(2); sc(i,2), y_lims(2); sc(i,2) y_lims(1)], ...
          'faces', [1, 2, 3, 4], ...
          'FaceColor', 'b', ...
          'FaceAlpha', 0.2,...
          'EdgeColor','none');
    th = text(midpt(sc(i,:)),y_lims(2)-0.125*(y_lims(2)-y_lims(1)),txt{i});
    th.HorizontalAlignment = 'center';
    th.VerticalAlignment = 'bottom';
end
hold off
ylabel('Decay corrected w/ Aus (Mafic)')
title('log_1_0(HP) vs. Age')

% Do the fft on the detrended data - zero pad to 20x length and then to
% nearest power of 2 (for speed)
pad_length = length(mean_hp_5) * 20;
pad_length = 2^(nextpow2(pad_length));

fft_5 = fft(mean_hp_5,pad_length);
Pyy = fft_5.*conj(fft_5)/pad_length;
f_5 = (1/dT)/pad_length*(0:floor(pad_length/2))';

% Find 3 highest peaks
[Peak, PeakIdx] = findpeaks(Pyy(1:floor(pad_length/2)+1));
[sorted_peak, sorted_peak_idx] = sort(Peak,'descend');
PeakIdx = PeakIdx(sorted_peak_idx);
sorted_peak = sorted_peak(1:3);
PeakIdx = PeakIdx(1:3);

% Plot the fft frequency curve in right plot
subplot(2,2,4)
plot(f_5,Pyy(1:floor(pad_length/2)+1))
for i = 1:length(PeakIdx)
    text(f_5(PeakIdx(i)), sorted_peak(i), sprintf('%6.2f Ma', 1./f_5(PeakIdx(i))));
end
title('Power spectral density')
xlabel('Frequency (1/Ma)')
xlim([0 0.05])




%% 6 Decay corrected, with Aus (Felsic)
figure()
% fft
mean_hp_6 = [];
mid_age = (agebins(2:end)+agebins(1:end-1))./2;


for i = 1:length(mid_age)
    % Get hp data - sio2 and decay corrected previously
    hp = data.hp_origin_nosio2(data.avg_age >= agebins(i) & data.avg_age < agebins(i+1) & data.sio2 > 60);
    % Convert to log scale
    hp = log10(hp);
    
    % Take mean
    mean_hp_6(i) = sum(hp)./length(hp);
end

% Find min and max age can interpolate/use for
findFirstLast = ~isnan(mean_hp_6);
lastIdx = find(findFirstLast(:),1,'last');
firstIdx = find(findFirstLast(:),1,'first');

% Modify mean_hp and mid_age to this
mid_age = mid_age(firstIdx:lastIdx);
mean_hp_6 = mean_hp_6(firstIdx:lastIdx);

% fft requires evenly spaced data - interp across NaN fields lost
mean_hp_6 = interp1(mid_age(~isnan(mean_hp_6)),mean_hp_6(~isnan(mean_hp_6)),mid_age);

% Remove mean/detrend
mean_hp_6 = detrend(mean_hp_6);
mean_hp_6 = detrend(mean_hp_6,0);

% Plot raw data in left subplot
subplot(2,2,1)
hold on
plot(mid_age,mean_hp_6,'-r')
hold off

y_lims = ylim;
hold on
sinfit = fit(mid_age,mean_hp_6,'sin1');
tx = text((max(mid_age) + min(mid_age))/2,y_lims(1)+0.125*(y_lims(2)-y_lims(1)),strcat('Period: ',num2str(1/(sinfit.b1/(2*pi)))));
plot(sinfit,'-k')
hold off
b = gca;
legend(b,'off');

hold on
sc = [173 335; 750 1071; 1500 1800; 2450 2720];
y_lims = ylim;
txt = {'Pangea', 'Rodina', 'Nuna', 'Kenorland'};
for i = 1:size(sc,1)
    p = patch('vertices', [sc(i,1), y_lims(1); sc(i,1), y_lims(2); sc(i,2), y_lims(2); sc(i,2) y_lims(1)], ...
          'faces', [1, 2, 3, 4], ...
          'FaceColor', 'b', ...
          'FaceAlpha', 0.2,...
          'EdgeColor','none');
    th = text(midpt(sc(i,:)),y_lims(2)-0.125*(y_lims(2)-y_lims(1)),txt{i});
    th.HorizontalAlignment = 'center';
    th.VerticalAlignment = 'bottom';
end
hold off
ylabel('SiO_2 + decay corrected w/ Aus (Felsic)')
title('log_1_0(HP) vs. Age')

% Do the fft on the detrended data - zero pad to 20x length and then to
% nearest power of 2 (for speed)
pad_length = length(mean_hp_6) * 20;
pad_length = 2^(nextpow2(pad_length));

fft_6 = fft(mean_hp_6,pad_length);
Pyy = fft_6.*conj(fft_6)/pad_length;
f_6 = (1/dT)/pad_length*(0:floor(pad_length/2))';

% Find 3 highest peaks
[Peak, PeakIdx] = findpeaks(Pyy(1:floor(pad_length/2)+1));
[sorted_peak, sorted_peak_idx] = sort(Peak,'descend');
PeakIdx = PeakIdx(sorted_peak_idx);
sorted_peak = sorted_peak(1:3);
PeakIdx = PeakIdx(1:3);

% Plot the fft frequency curve in right plot
subplot(2,2,2)
plot(f_6,Pyy(1:floor(pad_length/2)+1))
for i = 1:length(PeakIdx)
    text(f_6(PeakIdx(i)), sorted_peak(i), sprintf('%6.2f Ma', 1./f_6(PeakIdx(i))));
end
title('Power spectral density')
xlabel('Frequency (1/Ma)')
xlim([0 0.05])

%% 7 - Decay corrected, without Aus 1400-2000 Ma (Felsic)
mean_hp_7 = [];
mid_age = (agebins(2:end)+agebins(1:end-1))./2;
ind_aus = strcmpi(data.country,'AU') & data.avg_age < 2000 & data.avg_age > 1400;
ind_aus = ~ind_aus;


for i = 1:length(mid_age)
    % Get hp data - sio2 and decay corrected previously
    hp = data.hp_origin_nosio2(data.avg_age >= agebins(i) & data.avg_age < agebins(i+1) & ind_aus & data.sio2 > 60);
    % Convert to log scale
    hp = log10(hp);
    
    % Take mean
    mean_hp_7(i) = sum(hp)./length(hp);
end

% Find min and max age can interpolate/use for
findFirstLast = ~isnan(mean_hp_7);
lastIdx = find(findFirstLast(:),1,'last');
firstIdx = find(findFirstLast(:),1,'first');

% Modify mean_hp and mid_age to this
mid_age = mid_age(firstIdx:lastIdx);
mean_hp_7 = mean_hp_7(firstIdx:lastIdx);

% fft requires evenly spaced data - interp across NaN fields lost
mean_hp_7 = interp1(mid_age(~isnan(mean_hp_7)),mean_hp_7(~isnan(mean_hp_7)),mid_age);

% Remove mean/detrend
mean_hp_7 = detrend(mean_hp_7);
mean_hp_7 = detrend(mean_hp_7,0);

% Plot raw data in left subplot
subplot(2,2,3)
hold on
plot(mid_age,mean_hp_7,'-r')
hold off

y_lims = ylim;
hold on
sinfit = fit(mid_age,mean_hp_7,'sin1');
tx = text((max(mid_age) + min(mid_age))/2,y_lims(1)+0.125*(y_lims(2)-y_lims(1)),strcat('Period: ',num2str(1/(sinfit.b1/(2*pi)))));
plot(sinfit,'-k')
hold off
b = gca;
legend(b,'off');

hold on
sc = [173 335; 750 1071; 1500 1800; 2450 2720];
y_lims = ylim;
txt = {'Pangea', 'Rodina', 'Nuna', 'Kenorland'};
for i = 1:size(sc,1)
    p = patch('vertices', [sc(i,1), y_lims(1); sc(i,1), y_lims(2); sc(i,2), y_lims(2); sc(i,2) y_lims(1)], ...
          'faces', [1, 2, 3, 4], ...
          'FaceColor', 'b', ...
          'FaceAlpha', 0.2,...
          'EdgeColor','none');
    th = text(midpt(sc(i,:)),y_lims(2)-0.125*(y_lims(2)-y_lims(1)),txt{i});
    th.HorizontalAlignment = 'center';
    th.VerticalAlignment = 'bottom';
end
hold off
ylabel('Decay corrected w/ Aus (Felsic)')
title('log_1_0(HP) vs. Age')

% Do the fft on the detrended data - zero pad to 20x length and then to
% nearest power of 2 (for speed)
pad_length = length(mean_hp_7) * 20;
pad_length = 2^(nextpow2(pad_length));

fft_7 = fft(mean_hp_7,pad_length);
Pyy = fft_7.*conj(fft_7)/pad_length;
f_7 = (1/dT)/pad_length*(0:floor(pad_length/2))';

% Find 3 highest peaks
[Peak, PeakIdx] = findpeaks(Pyy(1:floor(pad_length/2)+1));
[sorted_peak, sorted_peak_idx] = sort(Peak,'descend');
PeakIdx = PeakIdx(sorted_peak_idx);
sorted_peak = sorted_peak(1:3);
PeakIdx = PeakIdx(1:3);

% Plot the fft frequency curve in right plot
subplot(2,2,4)
plot(f_7,Pyy(1:floor(pad_length/2)+1))
for i = 1:length(PeakIdx)
    text(f_7(PeakIdx(i)), sorted_peak(i), sprintf('%6.2f Ma', 1./f_7(PeakIdx(i))));
end
title('Power spectral density')
xlabel('Frequency (1/Ma)')
xlim([0 0.05])


%% 8 - Collision and subduction (Condie)
figure()
subplot(2,2,1)
[S,C,age] = orogen_hist(0,2200);

% Demean
S = detrend(S);
S = detrend(S,0);
plot(age,S,'-r')

y_lims = ylim;
hold on
sinfit = fit(age,S,'sin1');
tx = text((max(age) + min(age))/2,y_lims(1)+0.125*(y_lims(2)-y_lims(1)),strcat('Period: ',num2str(1/(sinfit.b1/(2*pi)))));
plot(sinfit,'-k')
hold off
b = gca;
legend(b,'off');


hold on
sc = [173 335; 750 1071; 1500 1800; 2450 2720];
y_lims = ylim;
txt = {'Pangea', 'Rodina', 'Nuna', 'Kenorland'};
for i = 1:size(sc,1)
    p = patch('vertices', [sc(i,1), y_lims(1); sc(i,1), y_lims(2); sc(i,2), y_lims(2); sc(i,2) y_lims(1)], ...
          'faces', [1, 2, 3, 4], ...
          'FaceColor', 'b', ...
          'FaceAlpha', 0.2,...
          'EdgeColor','none');
    th = text(midpt(sc(i,:)),y_lims(2)-0.125*(y_lims(2)-y_lims(1)),txt{i});
    th.HorizontalAlignment = 'center';
    th.VerticalAlignment = 'bottom';
end
hold off
ylabel('Subduction (Condie)')
title('Subduction vs age')

length_S = length(S)*20;
length_S = 2^nextpow2(length_S);

subplot(2,2,2)
S = fft(S,length_S);
Pyy = S.*conj(S)/length_S;
f_S = (1/dT)/length_S*(0:floor(length_S/2))';
[Peak, PeakIdx] = findpeaks(Pyy(1:floor(length_S/2)+1)); 
plot(f_S,Pyy(1:floor(length_S/2)+1))
% Get 3 highest peaks
[sorted_peak, sorted_peak_idx] = sort(Peak,'descend');
PeakIdx = PeakIdx(sorted_peak_idx);
sorted_peak = sorted_peak(1:3);
PeakIdx = PeakIdx(1:3);
for i = 1:length(PeakIdx)
    text(f_S(PeakIdx(i)), sorted_peak(i), sprintf('Peak = %6.6f Ma', 1./f_S(PeakIdx(i))));
end
%text(f(PeakIdx), Peak, sprintf('Peak = %6.3f', Peak));
title('Power spectral density')
xlabel('Frequency (Hz)')


subplot(2,2,3)

% Demean
C = detrend(C);
C = detrend(C,0);
plot(age,C,'-r')

y_lims = ylim;
hold on
sinfit = fit(age,C,'sin1');
tx = text((max(age) + min(age))/2,y_lims(1)+0.125*(y_lims(2)-y_lims(1)),strcat('Period: ',num2str(1/(sinfit.b1/(2*pi)))));
plot(sinfit,'-k')
hold off
b = gca;
legend(b,'off');



hold on
sc = [173 335; 750 1071; 1500 1800; 2450 2720];
y_lims = ylim;
txt = {'Pangea', 'Rodina', 'Nuna', 'Kenorland'};
for i = 1:size(sc,1)
    p = patch('vertices', [sc(i,1), y_lims(1); sc(i,1), y_lims(2); sc(i,2), y_lims(2); sc(i,2) y_lims(1)], ...
          'faces', [1, 2, 3, 4], ...
          'FaceColor', 'b', ...
          'FaceAlpha', 0.2,...
          'EdgeColor','none');
    th = text(midpt(sc(i,:)),y_lims(2)-0.125*(y_lims(2)-y_lims(1)),txt{i});
    th.HorizontalAlignment = 'center';
    th.VerticalAlignment = 'bottom';
end
hold off
ylabel('Collision (Condie)')
title('Collision vs age')

length_C = length(C)*20;
length_C = 2^nextpow2(length_C);

subplot(2,2,4)
C = fft(C,length_C);
Pyy = C.*conj(C)/length_C;
f_C = (1/dT)/length_C*(0:floor(length_C/2))';
[Peak, PeakIdx] = findpeaks(Pyy(1:floor(length_C/2)+1)); 
plot(f_C,Pyy(1:floor(length_C/2)+1))
% Get 3 highest peaks
[sorted_peak, sorted_peak_idx] = sort(Peak,'descend');
PeakIdx = PeakIdx(sorted_peak_idx);
sorted_peak = sorted_peak(1:3);
PeakIdx = PeakIdx(1:3);
for i = 1:length(PeakIdx)
    text(f_C(PeakIdx(i)), sorted_peak(i), sprintf('Peak = %6.6f Ma', 1./f_C(PeakIdx(i))));
end
%text(f(PeakIdx), Peak, sprintf('Peak = %6.3f', Peak));
title('Power spectral density')
xlabel('Frequency (Hz)')


%% 9 - eNd and eHf
figure()
agebins = (min_T:dT:max_T)';
filename = 'eNd_eHf.xlsx';
sheet = 1;

xlRange = 'C12:G10499';
Hfdata = xlsread(filename,sheet,xlRange);
xlRange = 'K13:M5154';
Nddata = xlsread(filename,sheet,xlRange);

medHf = [];
medNd = [];
Hfval = Hfdata(:,3);
Ndval = Nddata(:,2);
for i = 1:(length(agebins)-1)
    % Get age and hp data for age bin
    Hf = Hfval(Hfdata(:,1) >= agebins(i)-25 & Hfdata(:,1) < agebins(i+1)+25);
    Nd = Ndval(Nddata(:,1) >= agebins(i)-25 & Nddata(:,1) < agebins(i+1)+25);
    
    % Take median
    medHf(i) = median(Hf);
    medNd(i) = median(Nd);
end

xpoints_Hf = (agebins(2:end) + agebins(1:(end-1)))./2;
xpoints_Nd = (agebins(2:end) + agebins(1:(end-1)))./2;

% Find min and max age can interpolate/use for
findFirstLast = ~isnan(medHf);
lastIdx = find(findFirstLast(:),1,'last');
firstIdx = find(findFirstLast(:),1,'first');

% Modify mean_hp and mid_age to this
medHf = medHf(firstIdx:lastIdx);
xpoints_Hf = xpoints_Hf(firstIdx:lastIdx);

% Find min and max age can interpolate/use for
findFirstLast = ~isnan(medNd);
lastIdx = find(findFirstLast(:),1,'last');
firstIdx = find(findFirstLast(:),1,'first');

% Modify mean_hp and mid_age to this
medNd = medNd(firstIdx:lastIdx);
xpoints_Nd = xpoints_Nd(firstIdx:lastIdx);

% fft requires evenly spaced data - interp across NaN fields lost
medHf = interp1(xpoints_Hf(~isnan(medHf)),medHf(~isnan(medHf)),xpoints_Hf);
medNd = interp1(xpoints_Nd(~isnan(medNd)),medNd(~isnan(medNd)),xpoints_Nd);

% Remove mean/detrend
medHf = detrend(medHf);
medNd = detrend(medNd);

medHf = detrend(medHf,0);
medNd = detrend(medNd,0);

subplot(2,2,1)
plot(xpoints_Hf,medHf,'-r')

y_lims = ylim;
hold on
sinfit = fit(xpoints_Hf,medHf,'sin1');
tx = text((max(xpoints_Hf) + min(xpoints_Hf))/2,y_lims(1)+0.125*(y_lims(2)-y_lims(1)),strcat('Period: ',num2str(1/(sinfit.b1/(2*pi)))));
plot(sinfit,'-k')
hold off
b = gca;
legend(b,'off');



hold on
sc = [173 335; 750 1071; 1500 1800; 2450 2720];
y_lims = ylim;
txt = {'Pangea', 'Rodina', 'Nuna', 'Kenorland'};
for i = 1:size(sc,1)
    p = patch('vertices', [sc(i,1), y_lims(1); sc(i,1), y_lims(2); sc(i,2), y_lims(2); sc(i,2) y_lims(1)], ...
          'faces', [1, 2, 3, 4], ...
          'FaceColor', 'b', ...
          'FaceAlpha', 0.2,...
          'EdgeColor','none');
    th = text(midpt(sc(i,:)),y_lims(2)-0.125*(y_lims(2)-y_lims(1)),txt{i});
    th.HorizontalAlignment = 'center';
    th.VerticalAlignment = 'bottom';
end
hold off
ylabel('eHf')
title('eHf vs age')

subplot(2,2,3)
plot(xpoints_Nd,medNd,'-r')

y_lims = ylim;
hold on
sinfit = fit(xpoints_Nd,medNd,'sin1');
tx = text((max(xpoints_Nd) + min(xpoints_Nd))/2,y_lims(1)+0.125*(y_lims(2)-y_lims(1)),strcat('Period: ',num2str(1/(sinfit.b1/(2*pi)))));
plot(sinfit,'-k')
hold off
b = gca;
legend(b,'off');

hold on
sc = [173 335; 750 1071; 1500 1800; 2450 2720];
y_lims = ylim;
txt = {'Pangea', 'Rodina', 'Nuna', 'Kenorland'};
for i = 1:size(sc,1)
    p = patch('vertices', [sc(i,1), y_lims(1); sc(i,1), y_lims(2); sc(i,2), y_lims(2); sc(i,2) y_lims(1)], ...
          'faces', [1, 2, 3, 4], ...
          'FaceColor', 'b', ...
          'FaceAlpha', 0.2,...
          'EdgeColor','none');
    th = text(midpt(sc(i,:)),y_lims(2)-0.125*(y_lims(2)-y_lims(1)),txt{i});
    th.HorizontalAlignment = 'center';
    th.VerticalAlignment = 'bottom';
end
hold off
ylabel('eNd')
title('eNd vs age')

subplot(2,2,2)
length_Hf = length(medHf)*20;
length_Hf = 2^nextpow2(length_Hf);
Hf = fft(medHf,length_Hf);
Pyy = Hf.*conj(Hf)/length_Hf;
f_Hf = (1/dT)/length_Hf*(0:floor(length_Hf/2))';
[Peak, PeakIdx] = findpeaks(Pyy(1:floor(length_Hf/2)+1)); 
plot(f_Hf,Pyy(1:floor(length_Hf/2)+1))
% Get 3 highest peaks
[sorted_peak, sorted_peak_idx] = sort(Peak,'descend');
PeakIdx = PeakIdx(sorted_peak_idx);
sorted_peak = sorted_peak(1:3);
PeakIdx = PeakIdx(1:3);
for i = 1:length(PeakIdx)
    text(f_Hf(PeakIdx(i)), sorted_peak(i), sprintf('Peak = %6.6f Ma', 1./f_Hf(PeakIdx(i))));
end
%text(f(PeakIdx), Peak, sprintf('Peak = %6.3f', Peak));
title('Power spectral density')
xlabel('Frequency (Hz)')

subplot(2,2,4)
length_Nd = length(medNd)*20;
length_Nd = 2^nextpow2(length_Nd);
Nd = fft(medNd,length_Nd);
Pyy = Nd.*conj(Nd)/length_Nd;
f_Nd = (1/dT)/length_Nd*(0:floor(length_Nd/2))';
[Peak, PeakIdx] = findpeaks(Pyy(1:floor(length_Nd/2)+1)); 
plot(f_Nd,Pyy(1:floor(length_Nd/2)+1))
% Get 3 highest peaks
[sorted_peak, sorted_peak_idx] = sort(Peak,'descend');
PeakIdx = PeakIdx(sorted_peak_idx);
sorted_peak = sorted_peak(1:3);
PeakIdx = PeakIdx(1:3);
for i = 1:length(PeakIdx)
    text(f_Nd(PeakIdx(i)), sorted_peak(i), sprintf('Peak = %6.6f Ma', 1./f_Nd(PeakIdx(i))));
end
%text(f(PeakIdx), Peak, sprintf('Peak = %6.3f', Peak));
title('Power spectral density')
xlabel('Frequency (Hz)')


return
