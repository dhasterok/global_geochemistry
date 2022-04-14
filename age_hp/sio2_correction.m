function [hp_sio2,shift_var] = sio2_correction(sio2,hp,age,age_div,sio2_div,rock_origin,corr_val)
% hp_sio2 = sio2_correction(sio2,hp,age,age_div,sio2_div)
% Shifts HP to a common sio2 level - magnitude correction for rock type
% Weighted on number of samples per sio2 bin

fprintf('\n------------------\n')
fprintf('SiO_2 correction:\n')
fprintf('------------------\n\n')


%-------------------
% SiO2 fit
%-------------------

% SiO2 to shift all HP to
corr_val;

% Create output vector
hp_sio2 = nan(size(sio2));
shift_var = nan(size(sio2));

% r^2 value for raw data in log space:
r2 = corrcoef(sio2,log10(hp)).^2;
fprintf('r^2 value for raw HP to SiO2 in log space: %f\n',r2(1,2));


% Get median for each sio2 div
median_hp = [];
for i=1:(length(sio2_div)-1)
    ind = (sio2 >= sio2_div(i) & sio2 < sio2_div(i+1) & ~isnan(hp));
    median_hp(i) = median(log10(hp(ind)));
    median_sio2(i) = median(sio2(ind));
    no_data(i) = sum(ind);
end
weights_data = no_data ./ (sum(no_data));

r2 = corrcoef(median_sio2,median_hp).^2;
fprintf('r^2 value for median HP (every %d wt%%) to SiO2 in log space: %f\n',...
    sio2_div(2)-sio2_div(1),r2(1,2));

r2 = corrcoef(median_sio2(3:(end-1)),median_hp(3:(end-1))).^2;
fprintf('r^2 value for median HP (every %d wt%%) to SiO2 in log space (between 48 and 78): %f\n',...
    sio2_div(2)-sio2_div(1),r2(1,2));

% Create fit
x1 = median_sio2(:);
X = [ones(size(x1)) x1];
y = median_hp(:);
weights_data = weights_data(:);
newX = X;
for i = 1:size(X,2)
    newX(:,i) = weights_data.*X(:,i);
end
newY = y.*weights_data;
line_eq = (newX)\(newY);


% Estimate HP at sio2 designated using this equation
for i = 1:size(sio2,1)
    hp_sio2(i) = hp(i).*(10^((corr_val - sio2(i,:)).*line_eq(2)));
    shift_var(i) = (10^((corr_val - sio2(i,:)).*line_eq(2)));
end

% Plot data and equation
figure()
for i = 1:(length(sio2_div)-1)
    avg_sio2{i} = sio2(sio2 >= sio2_div(i) & sio2 < sio2_div(i+1) ...
        & ~isnan(hp) & ~isnan(sio2),:);
    avg_hp{i} = hp(sio2 >= sio2_div(i) & sio2 < sio2_div(i+1) ...
        & ~isnan(hp) & ~isnan(sio2),:);
end

sio2_div_2 = [sio2_div(1):((sio2_div(2)-sio2_div(1))*2):sio2_div(end)];
for i = 1:(length(sio2_div_2)-1)
    avg_sio2_2{i} = sio2(sio2 >= sio2_div_2(i) & sio2 < sio2_div_2(i+1) ...
        & ~isnan(hp) & ~isnan(sio2),:);
    avg_hp_2{i} = hp(sio2 >= sio2_div_2(i) & sio2 < sio2_div_2(i+1) ...
        & ~isnan(hp) & ~isnan(sio2),:);
end

[Qsio2,Qhp] = whisker(avg_sio2,avg_hp,'Color',[0.5 0.5 0.5],'Scale','log');
plot(0,0)

% Do stepped quantiles at half the resolution
[Qsio2_2,Qhp_2] = whisker(avg_sio2_2,avg_hp_2,'Color',[0.5 0.5 0.5],'Scale','log');
plot(0,0)

% Raw data
subplot(2,1,1)
% % Dont plot this temporarily
% hold on
% n = hist2d(sio2(~isnan(sio2) & ~isnan(hp)),...
%     log10(hp(~isnan(sio2) & ~isnan(hp))),...
%     sio2_div,-2:0.1:2);
% imagesc([min(sio2_div)+0.5 max(sio2_div)-0.5],[-2.5 2.5],n)
% colormap(flipud(bone));
% c = colorbar;
% c.Label.String = 'No. Samples';
% hold off
hpax([-2 2]);
xlim([min(sio2_div) max(sio2_div)])
caxis([0 400])


% Median points line is fit to (no weights shown)
hold on
%[Qsio2_2,Qhp_2] and sio2div_2

% Line 1 - 0.5 quantile
lx = length(Qsio2_2(:,3));
ly = size(Qhp_2,1);
X = zeros(2*lx,1);
Y = zeros(2*ly,5);
X(1:2:2*lx-1) = sio2_div_2(1:lx);
X(2:2:2*lx-2) = sio2_div_2(2:lx);
X(end)=sio2_div_2(end);
Y(1:2:2*ly-1,:) = Qhp_2(1:ly,:);
Y(2:2:2*ly,:)   = Qhp_2(1:ly,:);
plot(Qsio2_2(:,3),Qhp_2(:,3),'Color',[1 1 1],'Marker','.','MarkerSize',8,'linestyle', 'none');
median_plot = plot(X,Y(:,3),'linewidth',0.5,'Color',[1 1 1]);

% Line 2 - 0.05 quantile
median_plot = plot(X,Y(:,1),'linewidth',0.5,'Color',[0 0 0]);

% Line 3 - 0.25 quantile
median_plot = plot(X,Y(:,2),'linewidth',0.5,'Color',[0 0 0]);

% Line 4 - 0.75 quantile
median_plot = plot(X,Y(:,4),'linewidth',0.5,'Color',[0 0 0]);

% Line 5 - 0.95 quantile
median_plot = plot(X,Y(:,5),'linewidth',0.5,'Color',[0 0 0]);

hold off




% Best fit line
hold on
plot([40:5:85],line_eq(1)+line_eq(2).*[40:5:85],'-r')
hold off
golden
xlabel('SiO_2 [wt%]','FontSize',10)
ylabel('A [\muW m^{-3}]','FontSize',10)
set(gca,'Box','on');
fprintf('Line equation: log_1_0A = %.4f SiO_2 + %.4f\n',line_eq(2),line_eq(1))
xlim([40 85])

% RMSE of fit:
x_val = Qsio2(:,3);
y_val = Qhp(:,3);
y_val_est = line_eq(1)+line_eq(2).*Qsio2(:,3);
RMSE = 0;
for v = 1:length(Qsio2(:,3))
    RMSE = RMSE + weights_data(v).*((y_val_est(v) - y_val(v)).^2);
end
RMSE = sqrt(RMSE);
fprintf('Root mean square error (log units) on median points = %.2f\n',RMSE)


% Histogram of plutonic/volcanic
subplot(2,1,2)
histogram(sio2(strcmpi(rock_origin,'plutonic')),sio2_div,'DisplayStyle','stairs','EdgeColor','b');
hold on
histogram(sio2(strcmpi(rock_origin,'volcanic')),sio2_div,'DisplayStyle','stairs','EdgeColor','r');
histogram(sio2,sio2_div,'DisplayStyle','stairs','EdgeColor','k');
hold off
ylabel('No. data','FontSize',10);
xlabel('SiO2','FontSize',10);
xlim([40 85])
set(gca,'Box','on')
golden


return