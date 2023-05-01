function data = hpest_u_th(data,varargin)
% hpest_u_th - estimates heat production when U or Th are missing.
%
%   data = hpest_u_th(data);
%
%   Detailed plots may be added using the following command:
%   data = hpest_u_th(data,'Plot',true)

% parse inputs
p = inputParser;

addRequired(p,'data');
addParameter(p,'Plot',false,@islogical);

parse(p,data,varargin{:});

plotflag = p.Results.Plot;

% -----------------------------------
% all data irrespective of rock type
thind = data.th_ppm > 0;
uind = data.u_ppm > 0;

data.u_est = nan([height(data) 1]);
data.th_est = nan([height(data) 1]);

% compute ratio Th/U
data.ratio_th_u = nan([height(data) 1]);
data.ratio_th_u(uind & thind) = log10(data.th_ppm(uind & thind)./data.u_ppm(uind & thind));

% Treat igneous rocks separately from sedimentary rocks
% igneous rocks
% -----------------------------------
% For igneous rocks, use the median value to estimate U, Th
igind = rockgroup(data,'all igneous');
m_ig = quantile(data.ratio_th_u(uind & thind & igind),[0.025 0.25 0.50 0.75 0.975]);
sd_ig = std(data.ratio_th_u(uind & thind & igind));

% use Th to predict U
fprintf('Estimated %i U values (igneous)\n',sum(igind & thind & ~uind));
data.u_est(thind & igind) = data.th_ppm(thind & igind)/10^m_ig(3);

% use U to predict Th
fprintf('Estimated %i Th values (igneous)\n',sum(igind & uind & ~thind));
data.th_est(uind & igind) = 10^m_ig(3)*data.u_ppm(uind & igind);
fprintf('N: %i, Th/U: %.3f %.3f %.3f %.3f %.3f, mu: %f, sigma: %f\n\n',sum(uind & thind & igind),10.^m_ig,m_ig(3),sd_ig);


% sedimentary rocks
% -----------------------------------
% low sio2 --------------------------
% For sedimentary rocks, separate into carbonaceous and siliceous
sedind = rockgroup(data,'all seds');

% carbonate rocks are different than silicic rocks
sio2_low_ind = sedind & data.sio2 <= 30;
m_sio2_low = quantile(data.ratio_th_u(uind & thind & sio2_low_ind),[0.025 0.25 0.50 0.75 0.975]);
sd_sedlow = std(data.ratio_th_u(uind & thind & sio2_low_ind));
% use Th to predict U
fprintf('Estimated %i U values (sedimentary, SiO2 <= 30)\n',sum(sio2_low_ind & thind & ~uind));
data.u_est(thind & sio2_low_ind) = data.th_ppm(thind & sio2_low_ind)/10^m_sio2_low(3);

% use U to predict Th
fprintf('Estimated %i Th values (sedimentary, SiO2 <= 30)\n',sum(sio2_low_ind & uind & ~thind));
data.th_est(uind & sio2_low_ind) = 10^m_sio2_low(3)*data.u_ppm(uind & sio2_low_ind);
fprintf('N: %i, Th/U: %.3f %.3f %.3f %.3f %.3f, mu: %f, sigma: %f\n\n',sum(uind & thind & sio2_low_ind),10.^m_sio2_low,m_sio2_low(3),sd_sedlow);


% siliceous rocks
% high sio2 -------------------------
sio2_high_ind = sedind & data.sio2 >= 50;

m_sio2_high = quantile(data.ratio_th_u(uind & thind & sio2_high_ind),[0.025 0.25 0.50 0.75 0.975]);
sd_sedhigh = std(data.ratio_th_u(uind & thind & sio2_high_ind));
% use Th to predict U
fprintf('Estimated %i U values (sedimentary, SiO2 >= 50)\n',sum(sio2_high_ind & thind & ~uind));
data.u_est(thind & sio2_high_ind) = data.th_ppm(thind & sio2_high_ind)/10^m_sio2_high(3);

% use U to predict Th
fprintf('Estimated %i Th values (sedimentary, SiO2 >= 50)\n',sum(sio2_high_ind & uind & ~thind));
data.th_est(uind & sio2_high_ind) = 10^m_sio2_high(3)*data.u_ppm(uind & sio2_high_ind);
fprintf('N: %i, Th/U: %.3f %.3f %.3f %.3f %.3f, mu: %f, sigma: %f\n\n',sum(uind & thind & sio2_high_ind),10.^m_sio2_high,m_sio2_high(3),sd_sedhigh);


% mixed rocks
% mid sio2 --------------------------
sio2_mid_ind = sedind & 30 < data.sio2 & data.sio2 < 50;
m_sio2_mid = (m_sio2_high - m_sio2_low)/(50 - 30);
b_sio2_mid = m_sio2_low - 30*m_sio2_mid;

% use Th to predict U
fprintf('Estimated %i U values (sedimentary, 30 < SiO2 <= 50)\n',sum(sio2_mid_ind & thind & ~(data.u_ppm > 0)));
data.u_est(thind & sio2_mid_ind) = data.th_ppm(thind & sio2_mid_ind) ...
    ./ 10.^(m_sio2_mid(3)*data.sio2(thind & sio2_mid_ind) + b_sio2_mid(3));

% use U to predict Th
fprintf('Estimated %i Th values (sedimentary, 30 < SiO2 <= 50)\n',sum(sio2_mid_ind & uind & ~thind));
data.th_est(uind & sio2_mid_ind) = 10.^(m_sio2_mid(3)*data.sio2(uind & sio2_mid_ind) + b_sio2_mid(3)) ...
    .* data.u_ppm(uind & sio2_mid_ind);
fprintf('Th/U: %.3f %.3f %.3f %.3f %.3f\n\n',10.^m_sio2_mid);
fprintf('Th/U: %.3f %.3f %.3f %.3f %.3f\n\n',10.^b_sio2_mid);

% compute K2O from K
k2o = getK2O(data);

% predict heat production using estimated U, Th, when heat production is unknown
heat_production_mass = nan([height(data) 1]);
ind = k2o > 0 & thind & ~(uind);
heat_production_mass(ind) = computehp(k2o(ind),data.u_est(ind),data.th_ppm(ind));
ind = k2o > 0 & uind > 0 & ~(thind);
heat_production_mass(ind) = computehp(k2o(ind),data.u_ppm(ind),data.th_est(ind));

data.heat_production_est = heat_production_mass.*data.density_model;

% predict heat production using estimated U, Th, when heat production is known
ind = k2o > 0 & uind > 0 & thind > 0;
heat_production_u = nan([height(data) 1]);
heat_production_th = nan([height(data) 1]);
heat_production_u(ind) = data.density_model(ind).*computehp(k2o(ind),data.u_ppm(ind),data.th_est(ind));
heat_production_th(ind) = data.density_model(ind).*computehp(k2o(ind),data.u_est(ind),data.th_ppm(ind));


% make plots
% -----------------------------------
if plotflag
    figure;
    qrunavg(data.sio2(data.sio2 > 0 & uind & thind & igind), ...
        data.ratio_th_u(data.sio2 > 0 & uind & thind & igind),2);
    subplot(3,1,1:2);
    hpax([-2 2]);
    ylabel('Th/U Ratio');
    xlim([0 100]);
    title('Igneous & Metaigneous');
    for i = 1:5
        if i == 3
            plot([0 100],[m_ig(i) m_ig(i)],'r-','LineWidth',1.5);
        else
            plot([0 100],[m_ig(i) m_ig(i)],'g-','LineWidth',0.5);
        end
    end
    subplot(3,1,3);
    set(gca,'YScale','log');
    xlabel('SiO_2 (wt.%)');
    xlim([0 100]);


    figure;
    qrunavg(data.sio2(uind & thind & sedind), ...
        data.ratio_th_u(uind & thind & sedind),2);
    subplot(3,1,1:2);
    hpax([-2 2]);
    ylabel('Th/U Ratio');
    xlabel('SiO_2 (wt.%)');
    xlim([0 100]);
    title('Sedimentary & Metasedimentary');
    for i = 1:5
        if i == 3
            plot([0 30],[m_sio2_low(i) m_sio2_low(i)],'r-','LineWidth',1.5);
        else
            plot([0 30],[m_sio2_low(i) m_sio2_low(i)],'g-','LineWidth',0.5);
        end
    end
    for i = 1:5
        if i == 3
            plot([30 50],[m_sio2_mid(i)*30+b_sio2_mid(i) m_sio2_mid(i)*50+b_sio2_mid(i)],'r-','LineWidth',1.5);
        else
            plot([30 50],[m_sio2_mid(i)*30+b_sio2_mid(i) m_sio2_mid(i)*50+b_sio2_mid(i)],'g-','LineWidth',0.5);
        end
    end

    for i = 1:5
        if i == 3
            plot([50 100],[m_sio2_high(i) m_sio2_high(i)],'r-','LineWidth',1.5);
        else
            plot([50 100],[m_sio2_high(i) m_sio2_high(i)],'g-','LineWidth',0.5);
        end
    end
    subplot(3,1,3);
    set(gca,'YScale','log');
    ylabel('Th/U Ratio');
    xlim([0 100]);

    %return

    fprintf('Making igneous plot\n');
    figure_hpest(data(igind,:),m_ig,heat_production_u(igind),heat_production_th(igind));
    fprintf('Making sedimentary (SiO2 <= 50) plot\n');
    figure_hpest(data(sio2_low_ind,:),m_sio2_low,heat_production_u(sio2_low_ind),heat_production_th(sio2_low_ind));
    fprintf('Making sedimentary (SiO2 > 30) plot\n');
    figure_hpest(data(sio2_high_ind,:),m_sio2_high,heat_production_u(sio2_high_ind),heat_production_th(sio2_high_ind));
    fprintf('Making sedimentary (30 < SiO2 < 50) plot\n');
    figure_hpest(data(sio2_mid_ind,:),m_sio2_mid,heat_production_u(sio2_mid_ind),heat_production_th(sio2_mid_ind));
end

return



function figure_hpest(data,m,heat_production_u,heat_production_th)

uind = data.u_ppm > 0;
thind = data.th_ppm > 0;

fig1 = figure;

% uranium histograms
% --------------------
figure(fig1);
subplot(221); 
% observed uranium

histogram(log10(data.u_ppm(uind)),'BinWidth',0.1,'DisplayStyle','stairs');


hold on;
% modeled known uranium
histogram(log10(data.u_ppm(thind & uind)),'BinWidth',0.1,'DisplayStyle','stairs');
histogram(log10(data.u_est(thind & uind)),'BinWidth',0.1,'DisplayStyle','stairs');
% modeled unknown uranium
histogram(log10(data.u_est(thind & ~uind)),'BinWidth',0.1,'DisplayStyle','stairs');

legend('observed','known observed','known predicted','unknown predicted');
hpax([-3 3],'x');
xlabel('U (ppm)');
golden

fprintf('U misfit: %f\n',sqrt(sum(thind & uind)^-1*sum((log10(data.u_ppm(thind & uind)) - log10(data.u_est(thind & uind))).^2)));

% thorium histograms
% --------------------
figure(fig1);
subplot(223);
% observed thorium
histogram(log10(data.th_ppm(thind)),'BinWidth',0.1,'DisplayStyle','stairs');

hold on;
% modeled known thorium
histogram(log10(data.th_est(uind & thind)),'BinWidth',0.1,'DisplayStyle','stairs');
histogram(log10(data.th_est(uind & thind)),'BinWidth',0.1,'DisplayStyle','stairs');
% modeled unknown thorium
histogram(log10(data.th_est(uind & ~thind)),'BinWidth',0.1,'DisplayStyle','stairs');

legend('observed','known observed','known predicted','unknown predicted');
hpax([-3 3],'x');
xlabel('Th (ppm)');
golden

fprintf('Th misfit: %f\n',sqrt(sum(thind & uind)^-1*sum((log10(data.th_ppm(thind & uind)) - log10(data.th_est(thind & uind))).^2)));

% Th/U ratio
% --------------------
figure(fig1);
subplot(222);
histogram(data.ratio_th_u(uind),'BinWidth',0.1,'DisplayStyle','stairs');
hold on;
yl = get(gca,'YLim');
for i = 1:5
    if i == 3
        plot([m(i) m(i)],yl,'r-','LineWidth',1.5);
    else
        plot([m(i) m(i)],yl,'g-','LineWidth',0.5);
    end
end
hpax([-1 2],'x');
xlabel('Th/U ratio');
axis square;
golden

% 2-D histogram of Ratio Th/U
eu = [-3:0.05:3];
esi = [0:1:100];


fig2 = figure;
figure(fig2);
subplot(131);
er = [-3:0.05:3];
nbin = hist2d(data.sio2(uind),data.ratio_th_u(uind),esi,er);
imagesc(esi,er',log10(nbin));
colormap(parula);
caxis([-0.1 2]);
colorbar;

x = [0:5:100];
sio2 = data.sio2;
for i = 1:length(x)-1
    ind1 = (x(i) < sio2 & sio2 <= x(i+1));

    Rtemp{i} = data.ratio_th_u(ind1);
    Sitemp{i} = sio2(ind1);
end

[X,Y] = whisker(Sitemp,Rtemp,'Color',[0.5 0.5 0.5]);
%[X(:,3) 10.^Y(:,3)]

%hold on;
%p = polyfit(X(:,3),Y(:,3),5);
%plot(x,polyval(p,x),'-');
hold on;
for i = 1:5
    if i == 3
        plot([min(esi) max(esi)],[m(i) m(i)],'r-','LineWidth',1.5);
    else
        plot([min(esi) max(esi)],[m(i) m(i)],'g-','LineWidth',0.5);
    end
end
title(['Th/U ratio: Q2 (Q1,Q3) = ',num2str(10.^m(3)),' (',num2str(10.^m(2)),',',num2str(10.^m(4)),')']);

hpax([-1 2]);
ylabel('Th/U (ppm)');
xlabel('SiO_2 (wt.%)');
axis square;
axis xy;
set(gca,'Box','on');
%colorbar;


% heat production histograms
% --------------------
figure(fig1);
subplot(224);
% observed heat production

pos = data.heat_production > 0;
neg = data.heat_production < 0;

histogram(log10(data.heat_production(pos)),'BinWidth',0.1,'DisplayStyle','stairs');
hpax([-3 2],'x');
golden;


% --------------------
figure(fig2);
subplot(132);
nbin = hist2d(log10(data.u_ppm(uind)), ...
    log10(data.th_ppm(uind)) - m(3), ...
    eu,eu);
imagesc(eu,eu',nbin);
%colormap(flipud(gray));
colormap(parula);

hold on;
plot([-2 3],[-2 3],'-');
hpax([-2 3],'x');
xlabel('Observed U (ppm)');
hpax([-2 3],'y');
ylabel('Predicted U (ppm)');
axis square;
axis xy;
caxis([0 1000]);
%cbar;
colorbar;

subplot(133);
nbin = hist2d(log10(data.th_ppm(thind)), ...
    log10(data.u_ppm(thind)) - m(3), ...
    eu,eu);
imagesc(eu,eu',nbin);
%colormap(flipud(gray));
colormap(parula);

hold on;
plot([-2 3],[-2 3],'-');
hpax([-2 3],'x');
xlabel('Observed Th (ppm)');
hpax([-2 3],'y');
ylabel('Predicted Th (ppm)');
axis square;
axis xy;
caxis([0 1000]);
%cbar;
colorbar;


% heat production histogram
% --------------------
figure(fig1);
subplot(224);
hold on;
% modeled heat production
h = histogram(log10(heat_production_u),'BinWidth',0.1,'DisplayStyle','stairs');
h = histogram(log10(heat_production_th),'BinWidth',0.1,'DisplayStyle','stairs');
histogram(log10(data.heat_production_est(neg)),'BinWidth',0.1,'DisplayStyle','stairs');
histogram(log10([data.heat_production(pos);
    data.heat_production_est(data.heat_production_est > 0 & ~pos)]), ...
    'BinWidth',0.1,'DisplayStyle','stairs');
legend('observed','known predicted from U','known predicted Th','unknown predicted','observed + predicted');

return

