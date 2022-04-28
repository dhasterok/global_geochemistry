% PROCESSDB - A script to load and preprocess the global geochemical
%   database

% The script adds code subfolders that can be used to further analyze,
% process and plot the database.  It preprocesses the data using a specific
% set of options that allow for common treatment of all data.  Physical
% properties are then computed.

% start by clearing memory
%close all;  % if you want to get rid of figures too.
%clc;        % if you want to clear the command window
clear all;

% add directories with codes
addpath processing classification protolith physprop plotting toolbox  ...
    ternary fileio age_hp maptools worldgrid

% ----------------------------------------
% options and inputs
% ----------------------------------------

%str_input = input('Old tabledata.csv (y) or google drive (n)? y/n: ','s');

%if strcmpi(str_input,'y')
%    path = 'tabledata.csv';
%elseif strcmpi(str_input,'n')

%dbdate = '2018_06_20'; % granite paper
%dbdate = '2019_01_13'; % metamorphic protolith paper (contains antarctic
%data)
%dbdate = '2019_08_04'; % arc volcanics paper (japanese data, but not all antarctic data)
                        % bad country/oceanic locations
dbdate = '2022_02_03'; % update

% print a summary of properties by rock type at the end of this script?
output_stats = 1;

% note if you have downloaded this from GitHub, you may need to change this
% directory
%path = ['../database/export/',dbdate,'/database_',dbdate,'.csv'];
path = '/Users/dhasterok/Google Drive/My Drive/heat_production/database/export/';
       
% oxides for normalization
majors = {'sio2';'tio2';'al2o3';'feo_tot';'mgo';'cao';'na2o';'k2o';'p2o5'};
oxlist = {'SiO2'; 'TiO2'; 'Al2O3'; 'FeO'; 'MgO'; 'CaO'; 'Na2O'; 'K2O'; 'P2O5'};
tau = 0;                % oxide + volatile sum threshold for use
tau_vf = 0;             % oxide sum threshold for use

partial_flag = 1;       % Use partial chemical analyses? Yes = 1, No = 0

% REE summation scheme
% better to use reduced unless you are using very recent data (2018+) where
% all elements are routinely reported
%scheme = 'full'
scheme = 'reduced'; 

use_u_th = 0;           % Use (1) or exclude (0) h.p. estimates for missing
                        % U and Th

use_norm_oxide = 1;     % Pre-normalise (1) or not (0) major oxides

use_protolith_est = 1;  % Use (1) or exclude (0) estimated protoliths for
                        % metamorphic rocks with unknown protolith

% protolith classifier
%pclass = 'trained_WKNN_classifier_20190226';
pclass = 'trained_RUSBoost_Model_30l_1000s_20190226';

new_density_model = 0;  % Use (1) to update or (0) to use existing density
                        % models

% ----------------------------------------
% preparing data scripts
% ----------------------------------------
% read data file
fprintf('Reading database...\n');
data = readtable([path,dbdate,'/','database_',dbdate,'.csv']);

fprintf('Size of database: %i\n\n',height(data));

if strcmp(dbdate,'2018_06_20')
    data.rock_facies = data.rock_composition;
    data.rock_composition(:) = {''};
end

% reject minerals
data = data(~rockgroup(data,'mineral'),:);

%T = data;
% calculate total iron and iron ratio
fprintf('Convert all Fe to FeO and calculate Fe2+/Fe_total ratio...\n');
%sum(data{:,'fe2o3_tot'} > 0 | data{:,'fe2o3'} > 0 | data{:,'feo'} > 0 | data{:,'feo_tot'} > 0)
data = fefix(data);
%sum(data{:,'feo_tot'} > 0)


% convert cation data to oxides
fprintf('Convert cations to oxide data when missing...\n');
data = cat2ox(data);

% normalize oxide weights
if use_norm_oxide
    data = oxide_norm(data);
end

% remove partial analyses
if ~partial_flag
    ind = data.oxtotal > 80 & data.oxtotal < 120;
    data = data(ind,:);
end

if use_norm_oxide
    data = oxide_norm(data,oxlist);
end

% compute sum of REE's
fprintf('Computing sum of REEs from %s set...\n',scheme);
data = sumree(data,scheme);


% ----------------------------------------
% estimate ages
% ----------------------------------------
% in the future, this may be done to the native database so it won't need
% to be done here.
age_var = 200;
fprintf('Correcting ages...\n');
data.avg_age = age_correction(data,age_var);

%fprintf('Correcting oceanic and assigning oceanic ages...\n');
%data = correct_country(data);


% ----------------------------------------
% assign provinces
% ----------------------------------------
fprintf('Assigning provinces...\n');
[data,provinces] = addprovinces(data);

%figure;
%subplot(211);
%hold on;
%histogram(log10(provinces.area_km2));
%subplot(212);
%hold on;
%for i = 1:length(cont)
%    ind = strcmp(provinces.cont,cont{i});
%    histogram(log10(provinces.area_km2(ind)),'DisplayStyle','stairs');
%end
%legend(cont);
fprintf('\n');


% ----------------------------------------
% classifying rocks
% ----------------------------------------
%%%%% BEGIN TEMPORARY SECTION %%%%%
if isnumeric(data.rock_facies)
    data.rock_facies = [];
    data.rock_facies = data.rock_composition;
else
    ind = ~strcmp(data.rock_composition,'') & strcmp(data.rock_facies,'');
    data.rock_facies(ind) = data.rock_composition(ind);
end
data.rock_composition = [];
%%%%% END TEMPORARY SECTION %%%%%

% compute geochemical indicies and granite classes
fprintf('Computing geochemical indicies...\n');
data = geochem_index(data);

% compute O'Neill lambda values
fprintf('Computing O''Neill lambda values...\n');
data = lambdaree(data);

% protolith determination
% find metaigneous and metasedimentary rocks not labled as such based on
% rock_name
fprintf('Adjusting rock origins...\n');
data = adjust_origin(data);

if use_protolith_est
    fprintf('Classifying protolith using %s...\n',pclass);
    data.protolith_est = cell([height(data) 1]);
    data.protolith_est(:) = {''};

    [data.protolith_est,labels,protolithScore] = classifyProtolith(data,pclass);
else
    fprintf('Skipping protolith classification...\n',pclass);
end

% compute igneous names based on TAS (volcanic/plutonic) classification from Middlemost
% ESR 1994
fprintf('Computing TAS names...\n');
data.rock_type = cell([height(data) 1]);
data.rock_type(:) = {''};

ind = rockgroup(data,'igneous protolith');
data.rock_type(ind) = tas(data(ind,:),'Plot','none');

fprintf('Computing QAPF names...\n');
data.qap_name = cell([height(data) 1]);
data.qap_name(:) = {''};

%data.qap_name(ind) = QAP_plot(data(ind,:),'Norm','CIPW');

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

ind = rockgroup(data,'igneous protolith');
[data.sia_scheme(ind), ...
    data.frost_class1(ind),data.frost_class2(ind),data.frost_class3(ind)] ...
    = graniteclass(data(ind,:));

% compute sedimentary names from Mason (1967) Turekian (1968)
fprintf('Determining sedimentary classification...\n');
ind = rockgroup(data,'sedimentary protolith');
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
data(ind,:) = metamorphic_class(data(ind,:));


% ----------------------------------------
% computing physical properties
% ----------------------------------------
fprintf('\nComputing rock properties ...\n');
% estimate seismic velocities
fprintf('   P-velocity...');
data = vpest(data);
fprintf('done.\n');

fprintf('   S-velocity...');
data = vsest(data);
fprintf('done.\n');

% estimate densities
fprintf('   density...');
if new_density_model
    fprintf('updating models...');
    data = densest(data);
else
    fprintf('using existing models...');
    data = densest2(data);
end
fprintf('done.\n');

% estimate thermal conductivities
fprintf('   thermal_conductivity...');
data = tcest(data);
fprintf('done.\n');

% estimate heat production
fprintf('   heat production...');
data = hpest(data);
fprintf('done.\n\n');

% add heat production estimates for rocks without U or Th measurements
if use_u_th
    fprintf('Estimating heat production for missing Th, U...\n');
    data = hpest_u_th(data);
    ratio_vs_sio2(data);
    if use_u_th
        fprintf('\n   including estimated values in analysis...');
        ind = ~(data.heat_production > 0) & data.heat_production_est > 0;
        data.heat_production(ind) = data.heat_production_est(ind);
    else
        fprintf('\n   excluding estimated values in analysis...');
    end
    fprintf('done.\n\n');
else
    fprintf('Skipping Th/U heat production estimation...\n');
end

% estimate basalt liquidus
fprintf('Estimating basaltic liquidus...');
data.liquidus = nan([height(data) 1]);
ind = data.sio2 < 60 & rockgroup(data,'igneous protolith');
data.liquidus(ind) = basalt_liquidus(data(ind,:));
fprintf('done.\n\n');


% ----------------------------------------
% crustal thickness
% ----------------------------------------
%fprintf('Interpolate crustal thickness...\n');
%tmp = load('../data/seismic/crustal_thickness/crust1_crsthk.xyz');
%ct.longitude = unique(tmp(:,1));
%ct.latitude = unique(tmp(:,2));
%ct.moho = flipud(reshape(tmp(:,3),length(ct.longitude),length(ct.latitude))');
%            
%data.moho = nan([height(data),1]);
%
%data.moho = interp2(ct.longitude,ct.latitude,ct.moho,data.longitude,data.latitude);


% ----------------------------------------
% output stats
% ----------------------------------------
if output_stats
    data.loghp(data.heat_production > 0) = log(data.heat_production(data.heat_production > 0));
    rockstats = grpstats(data(data.heat_production > 0 & data.heat_production < 100,{'rock_type','sio2','density_model','p_velocity','s_velocity','thermal_conductivity','heat_production','loghp'}),'rock_type',{'mean','std'});
    writetable(rockstats,'rockstats.xlsx');
end
