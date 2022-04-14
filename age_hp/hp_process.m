close all;
clear all;
clc;

% ----------------------------------------
% options and inputs
% ----------------------------------------

use_u_th = 0; % Use (1) or exclude (0) h.p. estimates for missing U and Th
use_age_only = 0; % Use (1) or exclude (0)
use_dens_only = 0;
use_protolith_est = 0; % Use (1) or exclude (0) estimated protoliths for
    % metamorphic rocks with unknown protolith
    use_norm_oxide = 1; % Pre-normalise (1) or not (0) major oxides

% protolith classifier
pclass = 'protolithClassifier_2017_05_08_RUSboosted.mat';

%Load data
%dbdate = '2018_09_12';
%dbdate = '2018_06_20';
%dbdate = '2019_03_02';
%dbdate = '2018_09_12'; which one was i using..?
%dbdate = '2018_10_29';

dbdate = '2019_01_13'; % The data in the paper
% dbdate = '2018_09_12'; % The earlier version?

    path = ['../database/export/',dbdate,'/database_',dbdate,'.csv'];
%fmt = '%s%f%s%s%f%f%s%s%s%f%f%f%f%f%s%f%f%f%s%s%s%f%s%s%s%s%s%s%f%f%s%f%f%f%f%s%s%s%f%s%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f';
fmt = '%s%f%s%s%f%f%f%f%f%f%f%f%f%f%s%f%f%f%s%s%s%s%s%s%s%s%s%s%s%f%f%s%f%f%f%f%s%s%s%f%s%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f';


% dbdate = '2018_10_29';
% path = ['../database/export/',dbdate,'/database_',dbdate,'.csv'];
% fmt = ['%s%f%s%s%f%f%f%f%f%f%f%f%f%f%s%f%f%f%s%s%s%s%s%s%s%s%s%s%s%f%f',...
%     '%s%f%f%f%f%s%s%s%f%s%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f',...
%     '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f',...
%     '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f',...
%     '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f'];


    % read data file
fprintf('Reading database...\n');
if ~exist('data','var')
    data = readtable(path);
    %data = readtable(path,'Format',fmt);
end

if iscell(data.zr_ppm(1))
    zrind = cellfun(@isempty,data.zr_ppm);
    data.zr_ppm(zrind) = {'NaN'};
    data.zr_ppm = cellfun(@str2num,data.zr_ppm);
end

% Fix country
% Country
ind = strcmpi(data.country,'Algeria');
data.country(ind) = {'DZ'};
ind = strcmpi(data.country,'Antarctica');
data.country(ind) = {'AQ'};
ind = strcmpi(data.country,'Antarctica');
data.country(ind) = {'AQ'};
ind = strcmpi(data.country,'Argentina');
data.country(ind) = {'AR'};
ind = strcmpi(data.country,'Australia');
data.country(ind) = {'AU'};
ind = strcmpi(data.country,'Austria');
data.country(ind) = {'AT'};
ind = strcmpi(data.country,'brazil');
data.country(ind) = {'BR'};
ind = strcmpi(data.country,'Burkina Faso');
data.country(ind) = {'BF'};
ind = strcmpi(data.country,'Burundi');
data.country(ind) = {'BI'};
ind = strcmpi(data.country,'California');
data.country(ind) = {'US'};
ind = strcmpi(data.country,'Camaroon') | strcmpi(data.country,'Cameroon');
data.country(ind) = {'CM'};
ind = strcmpi(data.country,'Canada');
data.country(ind) = {'CA'};
ind = strcmpi(data.country,'Cape Verde');
data.country(ind) = {'CV'};
ind = strcmpi(data.country,'China');
data.country(ind) = {'CN'};
ind = strcmpi(data.country,'Congo');
data.country(ind) = {'CG'};
ind = strcmpi(data.country,'Czech Republic');
data.country(ind) = {'CZ'};
ind = strcmpi(data.country,'Dominican Republic');
data.country(ind) = {'DO'};
ind = strcmpi(data.country,'Ecuador');
data.country(ind) = {'EC'};
ind = strcmpi(data.country,'Egypt');
data.country(ind) = {'EG'};
ind = strcmpi(data.country,'England');
data.country(ind) = {'GB'};
ind = strcmpi(data.country,'Eritrea');
data.country(ind) = {'ER'};
ind = strcmpi(data.country,'Ethiopia');
data.country(ind) = {'ET'};
ind = strcmpi(data.country,'Finland');
data.country(ind) = {'FI'};
ind = strcmpi(data.country,'France');
data.country(ind) = {'FR'};
ind = strcmpi(data.country,'French Guyana');
data.country(ind) = {'GF'};
ind = strcmpi(data.country,'Germany');
data.country(ind) = {'DE'};
ind = strcmpi(data.country,'Greece');
data.country(ind) = {'GR'};
ind = strcmpi(data.country,'Greenland');
data.country(ind) = {'GL'};
ind = strcmpi(data.country,'Guinea');
data.country(ind) = {'GN'};
ind = strcmpi(data.country,'Hawaii');
data.country(ind) = {'US'};
ind = strcmpi(data.country,'India');
data.country(ind) = {'IN'};
ind = strcmpi(data.country,'Indonesia');
data.country(ind) = {'ID'};
ind = strcmpi(data.country,'Iran');
data.country(ind) = {'IR'};
ind = strcmpi(data.country,'Japan');
data.country(ind) = {'JP'};
ind = strcmpi(data.country,'Jordan');
data.country(ind) = {'JO'};
ind = strcmpi(data.country,'Kazakhstan');
data.country(ind) = {'KZ'};
ind = strcmpi(data.country,'Kyrgyzstan') | strcmpi(data.country,'Kyrgyzstan ');
data.country(ind) = {'KG'};
ind = strcmpi(data.country,'Madagascar');
data.country(ind) = {'MG'};
ind = strcmpi(data.country,'Mexico');
data.country(ind) = {'MX'};
ind = strcmpi(data.country,'Mongolia');
data.country(ind) = {'MN'};
ind = strcmpi(data.country,'Mozambique');
data.country(ind) = {'MZ'};
ind = strcmpi(data.country,'Namibia');
data.country(ind) = {'NA'};
ind = strcmpi(data.country,'Nigeria');
data.country(ind) = {'NG'};
ind = strcmpi(data.country,'Norway');
data.country(ind) = {'NO'};
ind = strcmpi(data.country,'Pakistan');
data.country(ind) = {'PK'};
ind = strcmpi(data.country,'Russia');
data.country(ind) = {'RU'};
ind = strcmpi(data.country,'Saudi Arabia');
data.country(ind) = {'SA'};
ind = strcmpi(data.country,'Scotland');
data.country(ind) = {'GB'};
ind = strcmpi(data.country,'Senegal');
data.country(ind) = {'SN'};
ind = strcmpi(data.country,'Sierra Leone');
data.country(ind) = {'SL'};
ind = strcmpi(data.country,'South Africa');
data.country(ind) = {'ZA'};
ind = strcmpi(data.country,'Spain');
data.country(ind) = {'ES'};
ind = strcmpi(data.country,'Sri Lanka');
data.country(ind) = {'LK'};
ind = strcmpi(data.country,'Sudan');
data.country(ind) = {'SD'};
ind = strcmpi(data.country,'Sweden');
data.country(ind) = {'SE'};
ind = strcmpi(data.country,'Switzerland');
data.country(ind) = {'CH'};
ind = strcmpi(data.country,'Tanzania');
data.country(ind) = {'TZ'};
ind = strcmpi(data.country,'Turkey');
data.country(ind) = {'TR'};
ind = strcmpi(data.country,'United States');
data.country(ind) = {'US'};
ind = strcmpi(data.country,'Vietnam');
data.country(ind) = {'VN'};
ind = strcmpi(data.country,'Yemen');
data.country(ind) = {'YE'};
ind = strcmpi(data.country,'Zambia');
data.country(ind) = {'ZM'};
ind = strcmpi(data.country,'Zimbabwe');
data.country(ind) = {'ZW'};
ind = strcmpi(data.country,'ocean');
data.country(ind) = {'Oceanic'};



















% reject minerals
data = data(~rockgroup(data,'mineral'),:);

% calculate total iron and iron ratio
fprintf('Convert all Fe to FeO and calculate Fe2+/Fe_total ratio...\n');
%data = fefix2(data);
data = fefix(data);

%Convert cation data to oxides
fprintf('Convert cations to oxide data when missing...\n');
data = cat2ox(data);
 
%Element correction: majors > 100, u/th > 1000, total > 110
data = element_correction(data);

%Normalize oxide weights
% oxides for normalization
majors = {'sio2';'tio2';'al2o3';'feo_tot';'mgo';'cao';'na2o';'k2o';'p2o5'};
oxlist = {'SiO2'; 'TiO2'; 'Al2O3'; 'Cr2O3'; 'FeO'; ...
    'MgO'; 'CaO'; 'Na2O'; 'K2O'; 'P2O5'; 'MnO'; 'NiO'; 'BaO'};
%oxlist = {'SiO2'; 'TiO2'; 'Al2O3'; 'FeO'; 'MgO'; 'CaO'; 'Na2O'; 'K2O'; 'P2O5'};
tau = 0; % oxide + volatile sum threshold for use
tau_vf = 0; % oxide sum threshold for use

if use_norm_oxide
    fprintf('Nomalizing to oxides: \n');
    for i = 1:length(oxlist)
        fprintf(' %s',oxlist{i});
    end
    fprintf('...\n');
    %data = oxide_norm(data,oxlist,tau,tau_vf);
    data = oxide_norm(data);
end


% ----------------------------------------
% classifying rocks
% ----------------------------------------
% compute geochemical indicies and granite classes
fprintf('Computing geochemical indicies...\n');
data = geochem_index(data);

% protolith determination
% find metaigneous and metasedimentary rocks not labled as such based on
% rock_name
fprintf('Adjusting rock origins...\n');
data = adjust_origin(data);

% fprintf('Classifying protolith using %s...\n',pclass);
% data.protolith_est = cell([height(data) 1]);
% data.protolith_est(:) = {''};
% 
% data.protolith_est = classifyProtolith(data,pclass);

% compute igneous names based on TAS (volcanic/plutonic) classification from Middlemost
% ESR 1994
fprintf('Computing TAS names...\n');
data.rock_type = cell([height(data) 1]);
data.rock_type(:) = {''};
ind = rockgroup(data,'all igneous');
%ind = rockgroup(data,'igneous protolith');

if use_age_only
    ind2 = (~isnan(data.age) | ~strcmpi(data.time_period,''));
    ind = and(ind,ind2);
end

if use_dens_only
    ind2 = (~isnan(data.density));
    ind = and(ind,ind2);
end

data.rock_type = cell([height(data) 1]);
data.rock_type(:) = {''};
data.rock_type(ind) = tas(data(ind,:),0);

fprintf('Computing QAPF names...\n');
data.qap_name = cell([height(data) 1]);
data.qap_name(:) = {''};
data.qap_name(ind) = QAP_plot(data(ind,:),'Norm','CIPW');

% compute SIA and Frost et al. (2001) granite classification schemes
fprintf('Determining igneous classification...\n');
data.sia_scheme = cell([height(data) 1]);
data.sia_scheme(:) = {''};
data.frost_class1 = cell([height(data) 1]);
data.frost_class2 = cell([height(data) 1]);
data.frost_class3 = cell([height(data) 1]);
data.frost_class1(:) = {''};
data.frost_class2(:) = {''};
data.frost_class3(:) = {''};
[data.sia_scheme(ind), ...
    data.frost_class1(ind),data.frost_class2(ind),data.frost_class3(ind)] ...
    = graniteclass(data(ind,:));

% compute sedimentary names from Mason (1967) Turekian (1968)
fprintf('Determining sedimentary classification...\n');
ind = rockgroup(data,'all seds');
data.quartz = nan([height(data) 1]);
data.feldspar = nan([height(data) 1]);
data.lithics = nan([height(data) 1]);
data(ind,:) = sedclass(data(ind,:));

% assign metamorphic facies and textures from rock names/descriptions
fprintf('Determining metamorphic facies...\n');
data.facies = cell([height(data) 1]);
data.facies(:) = {''};
data.texture = cell([height(data) 1]);
data.texture(:) = {''};

ind = rockgroup(data,'metamorphic');
if use_age_only
    ind2 = (~isnan(data.age) | ~strcmpi(data.time_period,''));
    ind = and(ind,ind2);
end

if use_dens_only
    ind2 = (~isnan(data.density));
    ind = and(ind,ind2);
end
data(ind,:) = metamorphic_class(data(ind,:));

%Adjust rock origin for metamorphic rocks based on name - look for granit
%etc.
data = adjust_origin(data);



% ----------------------------------------
% computing physical properties
% ----------------------------------------
fprintf('\nComputing rock properties ...\n');
% estimate seismic velocities
fprintf('   velocity...');
data = vpest(data);
fprintf('done.\n');


% estimate densities
fprintf('   density...');
data = densest(data);
fprintf('done.\n');


% estimate heat production
fprintf('   heat production...');
data = hpest(data);
fprintf('done.\n\n');


% add heat production estimates for rocks without U or Th measurements
% fprintf('Estimating heat production for missing U/Th...\n');
% data = hpest_u_th(data);
% ratio_vs_sio2(data);
% if use_u_th
%     fprintf('\n   including estimated values in analysis...');
%     ind = ~(data.heat_production > 0) & data.heat_production_est > 0;
%     data.heat_production(ind) = data.heat_production_est(ind);
% else
%     fprintf('\n   excluding estimated values in analysis...');
% end
% fprintf('done.\n\n');


% estimate basalt liquidus
fprintf('Estimating basaltic liquidus...');
data.liquidus = nan([height(data) 1]);
ind = data.sio2 < 60 & rockgroup(data,'all igneous');
data.liquidus(ind) = basalt_liquidus(data(ind,:));
fprintf('done.\n\n');