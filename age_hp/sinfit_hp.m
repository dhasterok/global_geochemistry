function sinfit_hp(age,hp,country,age_div)
% fft on hp vs age data to see if any dominant frequencies
% fft_hp(data2.avg_age,data2.hp_corrected,data2.country)

fprintf('\n----------------------\n')
fprintf('sinfit_hp: Best fit\n')
fprintf('----------------------\n\n')

% Proterozoic Australian data
ind_proto_aus = (strcmpi(country,'AU') & age >= 1400 & age <= 2000);

% Prepare binned plots:
for i = 1:length(age_div)-1
       ind = age_div(i) <= age & age < age_div(i+1);
       agebin.ind{i} = ind;
       avg_age_cell{i} = age(ind);
       hp_cell{i} = hp(ind);
       
       % Proto aus removed
       ind_noaus = (ind & ~ind_proto_aus);
       agebin2.ind{i} = ind_noaus;
       avg_age_cell_noaus{i} = age(ind_noaus);
       hp_cell_noaus{i} = hp(ind_noaus);
end

figure()
[agebin.Qage,agebin.Qhp] = whisker(avg_age_cell,hp_cell,'Color',[0.5 0.5 0.5],'Scale','log');
% No aus
[agebin2.Qage,agebin2.Qhp] = whisker(avg_age_cell_noaus,hp_cell_noaus,'Color',[0.5 0.5 0.5],'Scale','log');

% Not log
[agebin3.Qage,agebin3.Qhp] = whisker(avg_age_cell,hp_cell,'Color',[0.5 0.5 0.5]);
% Not log - no aus
[agebin4.Qage,agebin4.Qhp] = whisker(avg_age_cell_noaus,hp_cell_noaus,'Color',[0.5 0.5 0.5]);


% median_ages = agebin.Qage(:,3);
% median_hp = agebin.Qhp(:,3);
% % No aus
% median_ages_noaus = agebin2.Qage(:,3);
% median_hp_noaus = agebin2.Qhp(:,3);

% fprintf('Subtracted mean from data: Required for sin fit (no y shift in equation)\n')
% 
% % Detrend all hp
% agebin.Qhp(:,1) = agebin.Qhp(:,1) - nanmean(median_hp);
% agebin.Qhp(:,2) = agebin.Qhp(:,2) - nanmean(median_hp);
% agebin.Qhp(:,3) = agebin.Qhp(:,3) - nanmean(median_hp);
% agebin.Qhp(:,4) = agebin.Qhp(:,4) - nanmean(median_hp);
% agebin.Qhp(:,5) = agebin.Qhp(:,5) - nanmean(median_hp);
% % No aus
% agebin2.Qhp(:,1) = agebin2.Qhp(:,1) - nanmean(median_hp_noaus);
% agebin2.Qhp(:,2) = agebin2.Qhp(:,2) - nanmean(median_hp_noaus);
% agebin2.Qhp(:,3) = agebin2.Qhp(:,3) - nanmean(median_hp_noaus);
% agebin2.Qhp(:,4) = agebin2.Qhp(:,4) - nanmean(median_hp_noaus);
% agebin2.Qhp(:,5) = agebin2.Qhp(:,5) - nanmean(median_hp_noaus);
% 
% median_hp = agebin.Qhp(:,3) - nanmean(median_hp);
% % No aus
% median_hp_noaus = agebin2.Qhp(:,3) - nanmean(median_hp_noaus);



% HP at present day
subplot(1,2,1)
plot(0,0)
sqwavefill(agebin.Qhp,agebin.Qage(:,3),age_div,[0,0,0],'Test')
xlim([age_div(1) age_div(end)]);
set(gca,'Box','on');
axis square
hpax([-2 2]);

% fprintf('Best fit sin model:\n')
% [sinfit gof] = fit(median_ages,median_hp,'sin1')
% hold on
% plot(sinfit,'-r')
% hold off
% fprintf('Period: %.4f Ma\n',(2*pi)./sinfit.b1)
% fprintf('Amplitude (log-space): %.4f\n',sinfit.a1)
xlabel('Age [Ma]','FontSize',10)
ylabel('A [\muW m^{-3}]','FontSize',10)
title('Best fit sin wave')

% No aus
subplot(1,2,2)
plot(0,0)
sqwavefill(agebin2.Qhp,agebin2.Qage(:,3),age_div,[0,0,0],'Test')
xlim([age_div(1) age_div(end)]);
set(gca,'Box','on');
axis square
hpax([-2 2]);

% fprintf('Best fit sin model:\n')
% [sinfit_noaus gof_noaus] = fit(median_ages_noaus,median_hp_noaus,'sin1')
% hold on
% plot(sinfit_noaus,'-r')
% hold off
% fprintf('Period: %.4f Ma\n',(2*pi)./sinfit_noaus.b1)
% fprintf('Amplitude (log-space): %.4f\n',sinfit_noaus.a1)
xlabel('Age [Ma]','FontSize',10)
ylabel('A [\muW m^{-3}]','FontSize',10)
title('Best fit sin wave - no aus')




% HP at present day
subplot(1,2,1)
plot(0,0)
sqwavefill(agebin3.Qhp,agebin3.Qage(:,3),age_div,[0,0,0],'Test')
xlim([age_div(1) age_div(end)]);
set(gca,'Box','on');
axis square

% fprintf('Best fit sin model:\n')
% [sinfit gof] = fit(median_ages,median_hp,'sin1')
% hold on
% plot(sinfit,'-r')
% hold off
% fprintf('Period: %.4f Ma\n',(2*pi)./sinfit.b1)
% fprintf('Amplitude (log-space): %.4f\n',sinfit.a1)
xlabel('Age [Ma]','FontSize',10)
ylabel('A [\muW m^{-3}]','FontSize',10)

% No aus
subplot(1,2,2)
plot(0,0)
sqwavefill(agebin4.Qhp,agebin4.Qage(:,3),age_div,[0,0,0],'Test')
xlim([age_div(1) age_div(end)]);
set(gca,'Box','on');
axis square

% fprintf('Best fit sin model:\n')
% [sinfit_noaus gof_noaus] = fit(median_ages_noaus,median_hp_noaus,'sin1')
% hold on
% plot(sinfit_noaus,'-r')
% hold off
% fprintf('Period: %.4f Ma\n',(2*pi)./sinfit_noaus.b1)
% fprintf('Amplitude (log-space): %.4f\n',sinfit_noaus.a1)
xlabel('Age [Ma]','FontSize',10)
ylabel('A [\muW m^{-3}]','FontSize',10)



return
