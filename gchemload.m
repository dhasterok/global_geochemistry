function data = gchemload(varargin)
% GCHEMLOAD - Loads geochemical database
%
%   data = gchemload loads the most recent global whole-rock geochemical
%   database.
%
%   The database can be returned with the following option value pairs:
%
%       'Version'           select the database version (see below)
%       
%       'Normalization'     'none', 'anhydrous' (default) or 'hydrous',
%                           will normalize using a set of major oxides, 
%                           oxides can be set using 'Oxides'.  Note
%                           selecting hydrous will significantly reduce the
%                           size of the database.
%
%       'Oxides'            entered as a cell array of oxides.  The default
%                           list is {'SiO2'; 'TiO2'; 'Al2O3'; 'FeO'; 'MgO';
%                           'CaO'; 'Na2O'; 'K2O'; 'P2O5'}
%
%       'TotalTol'          tolerance for sum of oxides, i.e., sum of
%                           majors + loi/(h2o + co2 + so2 + f + cl) must be
%                           within 100 +/- tol in wt%.  Default is 10.
%
%       'REEScheme'         averaging scheme for rare earth elements,
%                           either 'full' or 'reduced' (default); see
%                           help sumree for details
%
%       'ProtolithClass'    true will predict the protolith of a
%                           metamorphic rock using a trained machine
%                           learning classifier (see Hasterok et al.
%                           (Computers & Geosci., 2019) for details;
%                           false (default)
%
%       'DerivedProperties' true (default) will compute a number of derived
%                           physical properties based on the chemistry.
%
%       'SpatialMetadata'   true (default) will add spatial metadata to the
%                           samples (e.g., crustal thickness)
%
%       'OutputStats'       true will output physical property statistics
%                           for each rock type; false (default)
%
%   Database versions:
%       '2022_02_03' (current)
%       '2019_08_04' (arc volcanics (japanese data, but not all antarctic
%                    data), bad country/oceanic locations
%       '2019_01_13' metamorphic protolith paper (contains antarctic data)
%       '2018_06_20' granite paper
%
% For additional information, contact:
%   Derrick Hasterok
%   Dept. Earth Sciences
%   Univ. of Adelaide
%   derrick.hasterok@adelaide.edu.au

addpath processing classification protolith physprop plotting toolbox  ...
    ternary fileio age_hp maptools worldgrid


% -----------------------------------------
% default options
% -----------------------------------------
% database version
dbdate = '2022_02_03';

% oxides for normalization
oxlist = {'SiO2'; 'TiO2'; 'Al2O3'; 'FeO'; 'MgO'; 'CaO'; 'Na2O'; 'K2O'; 'P2O5'};

% protolith classifier, see Hasterok et al. (Computers & Geosci., 2019)
%pclass = 'trained_WKNN_classifier_20190226';
pclass = 'trained_RUSBoost_Model_30l_1000s_20190226';

new_density_model = false; % update density models

tau = 10; % tolerance for oxide sum +/- volatiles (performed in oxide_norm.m)

% ----------------------------------------
% Inputs
% ----------------------------------------
p = inputParser;

addParameter(p,'Version',dbdate,@ischar);
addParameter(p,'Normalization','anhydrous',@ischar);
addParameter(p,'Oxides',oxlist,@iscell);
addParameter(p,'TotalTol',tau,@ischar);
addParameter(p,'REEScheme','reduced',@ischar);
addParameter(p,'ProtolithClass',true,@islogical);
addParameter(p,'OutputStats',false,@islogical);
addParameter(p,'DerivedProperties',true,@islogical);
addParameter(p,'SpatialMetadata',false,@islogical);

parse(p,varargin{:});

dbdate = p.Results.Version;
normtype = lower(p.Results.Normalization);
oxlist = p.Results.Oxides;
tau = p.Results.TotalTol;
scheme = p.Results.REEScheme;
use_protolith_est = p.Results.ProtolithClass;


% ----------------------------------------
% Load database
% ----------------------------------------

fmt = [];
% switch dbdate
%     case {'2017_05_08','2017_08_08','2017_09_22','2017_09_28','2017_10_13', ...
%             '2017_11_23','2018_04_27','2018_05_02','2018_06_08'}
%         fmt = '%s%f%s%s%f%f%s%s%s%f%f%f%f%f%s%f%f%f%s%s%s%f%s%s%s%s%s%s%f%f%s%f%f%f%s%s%s%f%s%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f';
%     case '2018_06_20'
%         fmt = '%s%f%s%s%f%f%s%s%s%f%f%f%f%f%s%f%f%f%s%s%s%f%s%s%s%s%s%s%f%f%s%f%f%f%f%s%s%s%f%s%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f';
%     case '2018_09_12'
%         fmt = '%s%f%s%s%f%f%f%f%f%f%f%f%f%f%s%f%f%f%s%s%s%s%s%s%s%s%s%s%s%f%f%s%f%f%f%f%s%s%s%f%s%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f';
%     case '2018_10_29'
%         fmt = '%s%f%s%s%f%f%f%f%f%f%f%f%f%f%s%f%f%f%s%s%s%s%s%s%s%s%s%s%s%f%f%s%f%f%f%f%s%s%s%f%s%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s%s%s%s%s';
%     case {'2019_01_13','2019_03_02'}
%         fmt = '%s%f%s%s%f%f%f%f%s%f%f%s%f%f%f%s%s%s%f%s%s%s%s%s%s%s%f%f%s%f%f%f%f%s%s%s%f%s%s%s%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s%s%s%s%s%s';
%     case {'2019_08_04'}
%         fmt = '%s%f%s%s%f%f%f%f%s%f%f%f%s%s%f%f%f%s%s%s%f%s%s%s%s%s%s%s%f%f%s%f%f%f%f%s%s%s%s%f%s%s%s%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s%s%s%s%s%s'
%     otherwise
%         fmt = [];
% end

path = ['../database/export/',dbdate,'/database_',dbdate,'.csv'];


% ----------------------------------------
% preparing data scripts
% ----------------------------------------
% read data file
fprintf('Reading database, %s...\n',dbdate);
if isempty(fmt)
    data = readtable(path);
else
    data = readtable(path,'Format',fmt);
end
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
switch normtype
    case 'none'
        fprintf('No normalization...\n');
    case 'anhydrous'
        fprintf('Normalizing to %s conditions...\n',normtype);
        data = oxide_norm(data,'Normalization',normtype,'Oxides',oxlist,'TotalTol',tau);
    case 'hydrous'
        fprintf('Normalizing to %s conditions...\n',normtype);
        % remove all data without LOI/volatiles (only accounting for H2O
        % and CO2)
        data = data(~isnan(data.loi) | ~any(~isnan(data{:,'h2o_plus','h2o','co2'})),:);
        data = oxide_norm(data,'Normalization',normtype,'Oxides',oxlist,'TotalTol',tau);
    otherwise
        error('Unknown norm type.');
end

% compute sum of REE's
if ~any(strcmp(scheme,{'full','reduced'}))
    error('Unknown REE summation scheme unknown.');
end
fprintf('Computing sum of REEs from %s set...\nSee sumree for details.\n',scheme);
data = sumree(data,scheme);


% ----------------------------------------
% estimate ages
% ----------------------------------------
age_var = 200;
fprintf('Correcting ages...\n');
data.avg_age = age_correction(data,age_var);


% ----------------------------------------
% assign provinces
% ----------------------------------------
% In the future this will be performed in SQL.
fprintf('Assigning provinces...\n');
[data,provinces] = addprovinces(data);
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

    [data.protolith_est,data.labels,data.protolithScore] = classifyProtolith(data,pclass);
else
    fprintf('Skipping protolith classification...\n');
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
if p.Results.DerivedProperties
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
    fprintf('Estimating heat production for missing Th, U...\n');
    data = hpest_u_th(data);
    ratio_vs_sio2(data);
    
    % estimate basalt liquidus
    fprintf('Estimating basaltic liquidus...');
    data.liquidus = nan([height(data) 1]);
    ind = data.sio2 < 60 & rockgroup(data,'igneous protolith');
    data.liquidus(ind) = basalt_liquidus(data(ind,:));
    fprintf('done.\n\n');
else
    fprintf('Skipping physical property estimates...\n');
end


% ----------------------------------------
% Spatial Metadata
% ----------------------------------------
if p.Results.SpatialMetadata
    % crustal thickness
    fprintf('Interpolate crustal thickness...\n');
    tmp = load('../data/seismic/crustal_thickness/crust1_crsthk.xyz');
    ct.longitude = unique(tmp(:,1));
    ct.latitude = unique(tmp(:,2));
    ct.moho = flipud(reshape(tmp(:,3),length(ct.longitude),length(ct.latitude))');

    data.moho = nan([height(data),1]);

    data.moho = interp2(ct.longitude,ct.latitude,ct.moho,data.longitude,data.latitude);
else
    fprintf('Skipping spatial metadata...\n');
end


% ----------------------------------------
% output stats
% ----------------------------------------
if p.Results.OutputStats
    fprintf('Writing rock statistics...\n');
    data.loghp(data.heat_production > 0) = log(data.heat_production(data.heat_production > 0));
    rockstats = grpstats(data(data.heat_production > 0 & data.heat_production < 100,{'rock_type','sio2','density_model','p_velocity','s_velocity','thermal_conductivity','heat_production','loghp'}),'rock_type',{'mean','std'});
    writetable(rockstats,'rockstats.xlsx');
end

fprintf('...Done!\n\n');

return