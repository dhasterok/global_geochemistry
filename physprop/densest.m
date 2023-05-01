function data = densest(data)
% DENSEST - estimates density using several methods
%
% Density is estimated using two methods:
%    (1) Christensen and Mooney (1995) where density is estimated from sesimic
%    velocity; and
%    (2) a model using a empirical density equation as a function of oxide
%    content.

% Behn & Kelemen use 7 oxides.  Normalize compositions to these 7 for
% density and velocity estimates.
oxides = {'SiO2'; 'TiO2'; 'Al2O3'; 'FeO'; 'MgO'; 'CaO'; 'Na2O'; 'K2O'; 'P2O5'};
tmp = oxide_norm(data,'Normalization','anhydrous','Oxides',oxides);
tmp = geochem_index(tmp);


% figure;
% subplot(221);
% plot(data.Fe_number,tmp.Fe_number,'.');
% subplot(222);
% plot(data.MALI,tmp.MALI,'.');
% subplot(223);
% plot(data.ASI,tmp.ASI,'.');
% subplot(224);
% plot(data.maficity,tmp.maficity,'.');

% estimation of density using Christensen and Mooney [JGR, 1995]
data.density_cm = 4929 - 13294./data.p_velocity;

% Note Behn and Kelemen (2003) did not produce a model in their paper, but did
% compute density as part of their velocity caluclations.


%rho = 60425 ...
%   - 600 * tmp.sio2 ...
%   - 374 * tmp.al2o3 ...
%   - 565 * tmp.mgo ...
%   - 474 * tmp.feo_tot ...
%   - 663 * tmp.cao ...
%   - 1498 * tmp.na2o ...
%   + 52 * tmp.k2o;

%rho = 3111 ...
%    - 4.576*tmp.sio2 ...
%    + 8.477*tmp.mgo ...
%    + 5.694*tmp.cao;

%rho = 2786.8 ...
%    - 1.5*tmp.sio2 ...
%    + 12.9*(tmp.mgo + tmp.feo_tot + tmp.cao) ...
%    - 15*tmp.k2o;

%rho = 4114 ...
%    - 15.4*tmp.sio2 ...
%    - 9.6*tmp.cao ...
%    - 67.2*tmp.na2o;

%rho = 4771 ...
%    - 23.2*tmp.sio2 ...
%    - 23.4*tmp.al2o3 ...
%    - 7.6*tmp.mgo ...
%    - 12.8*tmp.na2o;

%rho = 4751.3 ...
%    - 22.1*tmp.sio2 ...
%    - 8.7*tmp.mgo ...
%    - 18.7*tmp.cao ...
%    - 102.6*tmp.na2o;

%rho = 4215 ...
%    - 17*tmp.sio2 ...
%    - 14.7*tmp.al2o3 ...
%    - 1*tmp.cao ...
%    - 15.9*tmp.na2o;


%tmp = geochem_index(tmp);

%figure;
%subplot(221);
%plot(data.Fe_number,tmp.Fe_number,'.');
%title('Fe_number');
%xl = get(gca,'XLim');
%hold on;
%plot(xl,xl,'-');
%
%subplot(222);
%plot(data.MALI,tmp.MALI,'.');
%title('MALI');
%xl = get(gca,'XLim');
%hold on;
%plot(xl,xl,'-');
%
%subplot(223);
%plot(data.ASI,tmp.ASI,'.');
%title('ASI');
%xl = get(gca,'XLim');
%hold on;
%plot(xl,xl,'-');
%
%subplot(224);
%plot(data.maficity,tmp.maficity,'.');
%title('maficity');
%xl = get(gca,'XLim');
%hold on;
%plot(xl,xl,'-');
%
% without garnet geochemical index model
% A shift may be applied to fit density observations.  For example comparison to
% PETROCH specific gravity measurements, the shift is -40 kg m^-3 for plutonic
% rocks and -51 kg m^-3 for volcanic rocks.  A shift of -45 kg m^-3 is quite
% reasonable for both.
%shift = -45
shift = -120;
rho = shift + 2606.5 + 174.7*tmp.Fe_number - 12.0*tmp.MALI + 49.6*tmp.ASI + 636.0*tmp.maficity;

% only use data with reasonable oxide percentages
ind = (tmp.sio2 > 100 | ...
    tmp.tio2 > 100 | ...
    tmp.al2o3 > 100 | ...
    tmp.feo_tot > 100 | ...
    tmp.mgo > 100 | ...
    tmp.cao > 100 | ...
    tmp.na2o > 100 | ...
    tmp.k2o > 100 | ...
    tmp.p2o5 > 100);

rho(ind) = NaN;

% remove outliers
rho(rho<2400 | rho>3600) = NaN;
rho = rho(:);

data.density_bk = rho;

fprintf('\n  BK Density Model\n');
fprintf('  ------------------\n');
density_obs = data.density;
density_obs(density_obs < 10) = density_obs(density_obs < 10)*1000;
data.density = density_obs;

ind = ~strcmp('carbonatite',data.rock_type) & rockgroup(data,'igneous') & data.mgo < 18;

ax = [2400 3600];


% indicies into different groups for density analysis
igind = rockgroup(data,'igneous protolith');
sedind = rockgroup(data,'sedimentary protolith');

% Group 3 - Igneous carbonatites
ind3 = strcmp('carbonatite',data.rock_type);
% Group 4 - Sedimentary carbonates
ind4 = strcmp('dolomite',data.rock_type) | strcmp('limestone',data.rock_type);
% Group 1 - Low-magnesian, igneous and sedimentary rocks
%   sans carbonates/carbonatites
ind1 = ((data.mgo < 18 & rockgroup(data,'all igneous')) | rockgroup(data,'all seds')) & ~(ind3 | ind4);
% Group 2 - High-magnesian, igneous and sedimentary rocks
%   sans carbonates/carbonatites
ind2 = (data.mgo >= 18 & rockgroup(data,'all igneous')) & ~(ind3 | ind4);

figure;
msize = 3;
% group 1
plot(data.sio2(ind1 & igind),density_obs(ind1 & igind),'o','MarkerSize',msize);
hold on;
plot(data.sio2(ind1 & sedind),density_obs(ind1 & sedind),'d','MarkerSize',msize);
% group 2
plot(data.sio2(ind2 & igind),density_obs(ind2 & igind),'o','MarkerSize',msize);
plot(data.sio2(ind2 & sedind),density_obs(ind2 & sedind),'d','MarkerSize',msize);
% group 4
plot(data.sio2(ind4),density_obs(ind4),'d','MarkerSize',msize);
% group 3
plot(data.sio2(ind3),density_obs(ind3),'o','MarkerSize',msize);
xlabel('SiO_2 (wt.%)');
ylabel('Observed Density (wt.%)');
ylim([2000 4000])
legend('Group 1 - Low-MgO volcanic & plutonic', ...
    'Group 1 - sedimentary (non-carbonate)', ...
    'Group 2 - High-MgO volcanic', ...
    'Group 4 - Carbonates', ...
    'Group 3 - Carbonatite');

% figure;
% msize = 2;
% el = {'tio2','al2o3','feo_tot','mgo','cao','mno','na2o','k2o','p2o5', ...
% 'Mg_number','Fe_number','MALI','ASI','maficity','CIA','WIP','spar', ...
% 'qtzindex','R1','R2','CPA'};
% for i = 1:length(el)
%     if i > 9
%         if i == 10
%             figure;
%         end
%         subplot(3,4,i-9);
%     else
%         subplot(3,3,i);
%     end
%     plot(data{ind1 & igind,el{i}},density_obs(ind1 & igind),'o','MarkerSize',msize);
%     hold on;
%     plot(data{ind1 & sedind,el{i}},density_obs(ind1 & sedind),'d','MarkerSize',msize);
%     % group 2
%     plot(data{ind2 & igind,el{i}},density_obs(ind2 & igind),'o','MarkerSize',msize);
%     plot(data{ind2 & sedind,el{i}},density_obs(ind2 & sedind),'d','MarkerSize',msize);
%     % group 4
%     plot(data{ind4,el{i}},density_obs(ind4),'d','MarkerSize',msize);
%     % group 3
%     plot(data{ind3,el{i}},density_obs(ind3),'o','MarkerSize',msize);
%     xlabel([el{i},' (wt.%)']);
%     ylabel('Observed Density (wt.%)');
%     ylim([2000 4000])
% end

ind = ind & ~isnan(data.density) & ~isnan(data.density_bk);
misfit = sqrt(sum((density_obs(ind) - data.density_bk(ind)).^2)/sum(ind));
fprintf('  Misfit: %.0f\n',misfit);

figure;
subplot(121); hold on;
plot(density_obs(ind),data.density_bk(ind),'.');
plot(ax,ax,'-');
plot(ax,ax+100,'--',ax,ax-100,'--');
axis([ax ax]);
axis square;
ylabel('BK Density Model [kg m^{-3}]');
xlabel('Observed Density [kg m^{-3}]');
set(gca,'Box','on');

subplot(122);
histogram(data.density_bk(ind) - density_obs(ind),'BinWidth',20);
xlim([-400 400]);
xlabel('Density Error (kg m^{-3}');
axis square;



% density computed from SiO2
% ---------------------------
%terms = {'sio2','sio2^-1','mgo','cao','al2o3','na2o'};
%x = [ones(size(data.sio2)) data.sio2 data.sio2.^-1 data.mgo data.cao data.al2o3 data.na2o];
terms = {'Fe_number','maficity','MALI','loi'};
x = [ones(size(data.sio2)) data.Fe_number data.maficity data.MALI data.loi];
[data.density_sio2_loi,m_sio2_loi,mint_sio2_loi] = compute_density(ind1 & data.sio2 > 0,data,density_obs,x,terms,'SiO2 w/LOI');
data.density_sio2_loi = x*m_sio2_loi;

terms = {'Fe_number','maficity','MALI'};
x = [ones(size(data.sio2)) data.Fe_number data.maficity data.MALI];
[data.density_sio2,m_sio2,mint_sio2] = compute_density(ind1 & data.sio2 > 0,data,density_obs,x,terms,'SiO2');
data.density_sio2 = x*m_sio2;


% density of high-Mg igneous rocks
% ---------------------------
%terms = {'Mg_number','sio2','maficity','cao','al2o3'};
terms = {'mgo','cao','loi'};
[data.density_himg_loi,m_himg_loi,mint_himg_loi] = compute_density(ind2,data,density_obs,[],terms,'High-Mg w/LOI');

terms = {'mgo','cao'};
[data.density_himg,m_himg,mint_himg] = compute_density(ind2,data,density_obs,[],terms,'High-Mg');

% density of carbonatites
% ---------------------------
terms = {'maficity','feo_tot','p2o5','loi'};
[data.density_igcarb_loi,m_igcarb_loi,mint_igcarb_loi] = compute_density(ind3,data,density_obs,[],terms,'Carbonatite w/LOI');

terms = {'feo_tot','p2o5'};
[data.density_igcarb,m_igcarb,mint_igcarb] = compute_density(ind3,data,density_obs,[],terms,'Carbonatite');


% density of carbonates
% ---------------------------
terms = {'sio2','cao','mgo','loi'};
[data.density_sedcarb_loi,m_sedcarb_loi,mint_sedcarb_loi] = compute_density(ind4,data,density_obs,[],terms,'Carbonate w/LOI');

terms = {'sio2','cao','mgo'};
[data.density_sedcarb,m_sedcarb,mint_sedcarb] = compute_density(ind4,data,density_obs,[],terms,'Carbonate');

data.density_model = data.density_bk;
%data.density_model(ind1 & isnan(data.density_bk)) = data.density_sio2(ind1 & isnan(data.density_bk)); 
data.density_model(ind1) = data.density_sio2(ind1);
data.density_model(ind2) = data.density_sio2(ind2);
%data.density_model(ind2) = data.density_himg(ind2);
data.density_model(ind3) = data.density_igcarb(ind3);
data.density_model(ind4) = data.density_sedcarb(ind4);

data.density_model(data.density_model < 2400 | data.density_model > 3600) = NaN;

figure;
msize = 3;
% group 1
plot(data.sio2(ind1 & igind & density_obs > 0), ...
    data.density_model(ind1 & igind & density_obs > 0), ...
    'o','MarkerSize',msize);
hold on;
plot(data.sio2(ind1 & sedind & density_obs > 0), ...
    data.density_model(ind1 & sedind & density_obs > 0), ...
    'd','MarkerSize',msize);

% group 2
plot(data.sio2(ind2 & igind & density_obs > 0), ...
    data.density_model(ind2 & igind & density_obs > 0), ...
    'o','MarkerSize',msize);
plot(data.sio2(ind2 & sedind & density_obs > 0), ...
    data.density_model(ind2 & sedind & density_obs > 0), ...
    'd','MarkerSize',msize);

% group 4
plot(data.sio2(ind4 & density_obs > 0), ...
    data.density_model(ind4 & density_obs > 0), ...
    'd','MarkerSize',msize);

% group 3
plot(data.sio2(ind3 & density_obs > 0), ...
    data.density_model(ind3 & density_obs > 0), ...
    'o','MarkerSize',msize);

xlabel('SiO_2 (wt.%)');
ylabel('Modeled Density (wt.%)');
ylim([2000 4000])
legend('Group 1 - Low-MgO volcanic & plutonic', ...
    'Group 1 - sedimentary (non-carbonate)', ...
    'Group 2 - High-MgO volcanic', ...
    'Group 4 - Carbonates', ...
    'Group 3 - Carbonatite');


figure;
plot(data.density(ind1 & igind),data.density_model(ind1 & igind), ...
    'o','MarkerSize',msize);
hold on;
plot(data.density(ind1 & sedind),data.density_model(ind1 & sedind), ...
    'd','MarkerSize',msize);
plot(data.density(ind2 & igind),data.density_model(ind2 & igind), ...
    'o','MarkerSize',msize);
plot(data.density(ind2 & sedind),data.density_model(ind2 & sedind), ...
    'd','MarkerSize',msize);
plot(data.density(ind4),data.density_model(ind4), ...
    'd','MarkerSize',msize);
plot(data.density(ind3),data.density_model(ind3), ...
    'o','MarkerSize',msize);
ax = [2400 3600];
plot(ax,ax,'k-');
axis([ax ax]);
xlabel('Observed Density [kg m^{-3}]');
ylabel('Modeled Density [kg m^{-3}]');
legend('Group 1 - Low-MgO volcanic & plutonic', ...
    'Group 1 - sedimentary (non-carbonate)', ...
    'Group 2 - High-MgO volcanic', ...
    'Group 4 - Carbonates', ...
    'Group 3 - Carbonatite', ...
    'Location','EastOutside');
axis square;

%if saveflag
%    addpath latex;
%    
%    fid = fopen('density_models.tex','w');
%    
%    tableheader(fid, ...
%        'Density model.', ...
%        'summary', ...
%        'lccccc');
%    fprintf(fid,'        {} & Low-Mg & Low-Mg & High-Mg & High-Mg & igneous & igneous & sedimentary & sedimentary \\\\\n');
%    fprintf(fid,'        Parameter & silicates & silicates & silicates & silicates & carbonatite & carbonatite & carbonate & carbonate \\\\\n');
%    hline(fid);
%    fprintf(fid,'        %.1f (%.1f / %.1f) & ',m_sio2(1),mint_sio2(1,1),mint_sio2(1,2));
%    fprintf(fid,'%.1f (%.1f / %.1f) & ',m_sio2_loi(1),mint_sio2_loi(1,1),mint_sio2_loi(1,2));
%    fprintf(fid,'%.1f (%.1f / %.1f) & ',m_himg(1),mint_himg(1,1),mint_himg(1,2));
%    fprintf(fid,'%.1f (%.1f / %.1f) & ',m_himg_loi(1),mint_himg_loi(1,1),mint_himg_loi(1,2));
%    fprintf(fid,'%.1f (%.1f / %.1f) & ',m_igcarb(1),mint_igcarb(1,1),mint_igcarb(1,2));
%    fprintf(fid,'%.1f (%.1f / %.1f) & ',m_igcarb_loi(1),mint_igcarb_loi(1,1),mint_igcarb_loi(1,2));
%    fprintf(fid,'%.1f (%.1f / %.1f) & ',m_sedcarb(1),mint_sedcarb(1,1),mint_sedcarb(1,2));
%    fprintf(fid,'%.1f (%.1f / %.1f) & ',m_sedcarb(1),mint_sedcarb(1,1),mint_sedcarb(1,2));
%    fprintf(fid,'        misfit & \\\\\n');
%    hline(fid);
%    
%    
%    hline(fid);
%    tablefooter(fid);
%    fprintf('\n\n');
%end

return


function [density_model,m,mint] = compute_density(ind,data,density_obs,x,terms,title)

ind_fit = ind & ~isnan(density_obs);
density_model = nan([height(data) 1]);

if isempty(x)
    x = [ones([height(data) 1]) data{:,terms}];
end

[m,mint,r,rint,stats] = regress(density_obs(ind_fit),x(ind_fit,:));

fprintf('\n  %s Density Model\n',title);
fprintf('  -------------------------\n');
fprintf('  rho_0:   %f (%f,%f)\n',m(1),mint(1,:));
for i = 1:length(terms)
    fprintf('  %-7s:   %f (%f,%f)\n',terms{i},m(i+1),mint(i+1,:));
end

%eig(x(2:end,:)'*x(2:end,:))

density_model(ind) = x(ind,:)*m;

ind_fit = ind_fit & ~isnan(density_model);
misfit = sqrt( sum( (density_obs(ind_fit) - density_model(ind_fit)).^2 ) ...
    / sum(ind_fit) );
fprintf('  N:      %i\n',sum(ind_fit));
fprintf('  Misfit: %.0f\n',misfit);

figure;
ax = [2400 3600];

subplot(121); hold on;
plot(density_obs(ind_fit),density_model(ind_fit),'.');
plot(ax,ax,'-');
plot(ax,ax+100,'--',ax,ax-100,'--');

axis([ax ax]);
axis square;
ylabel([title,' Density Model [kg m^{-3}]']);
xlabel('Observed Density [kg m^{-3}]');
set(gca,'Box','on');

subplot(122);
histogram(density_model(ind_fit) - density_obs(ind_fit),'BinWidth',20);
xlim([-400 400]);
xlabel('Density Error [kg m^{-3}]');
axis square;

return
