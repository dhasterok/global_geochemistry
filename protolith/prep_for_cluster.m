function [traindata,testdata,scores] = prep_for_cluster(data,varargin);
% PREP_FOR_CLUSTER - Creates training and test datasets for cluster analysis.
%
%   [traindata,testdata] = prep_for_cluster(data) produces a training
%   dataset (~90% of data) and testing dataset (~10% of data) to train a
%   classifier for metamorphic protoliths.
%
%   [traindata,testdata] = prep_for_cluster(data,option) allows for chaning
%   data used for training.  option == 1 uses only igneous and sedimentary
%   rocks, option == 2 uses only metasedimentary and metaigneous rocks.
%   All other options use the default (metaigneous and igneous vs
%   metasedimentary and sedimentary).

try addpath coda
catch
    try addpath ../coda
    catch error('Could not find CoDA directory.');
    end
end    

% fields used for training
major = {'sio2','tio2','al2o3','feo_tot','mgo','cao','na2o','k2o','p2o5'};
trace = {'cr2o3','mno','nio','bao','sro','cs_ppm','rb_ppm','co_ppm', ...
    'zn_ppm','cu_ppm','v_ppm','ga_ppm','pb_ppm','sc_ppm','y_ppm', ...
    'th_ppm','u_ppm','zr_ppm','hf_ppm','nb_ppm','ta_ppm'};
ree = {'la_ppm','ce_ppm','pr_ppm','nd_ppm','sm_ppm','eu_ppm','gd_ppm', ...
    'tb_ppm','dy_ppm','ho_ppm','er_ppm','tm_ppm','yb_ppm','lu_ppm'};
%index = {'Fe_number','MALI','ASI','maficity','CIA','spar','qtzindex','CAI','AI','CPA'};

% determine option
% default values
option = 3; % input data
equal = 0; % flag for partitioning database
reserve = 0.1; % reserve 10% of the data for a postiori
usenan = 0;
fieldsets = {'major'};
centering_fcn = 'none';
pcaflag = 0;

opt = 1;
if nargin > 1
    while opt < nargin
        switch lower(varargin{opt})
            case 'subset'
                option = varargin{opt+1};
                opt = opt + 2;
            case 'equal' % default unequal (0)
                equal = 1;
                opt = opt + 1;
            case 'reserve' % 10%
                reserve = varargin{opt+1};
                opt = opt + 2;
            case 'fields'
                fieldsets = varargin{opt+1};
                opt = opt + 2;
            case 'usenan'
                usenan = 1;
                opt = opt + 1;
            case 'pca'
                pcaflag = 1;
                opt = opt + 1;
            case 'recenter'
                centering_fcn = varargin{opt+1};
                opt = opt + 2;
                if ~any(strcmp(centering_fcn,{'none','clr','ilr'}))
                    error('Unknown centering function.');
                end
            otherwise
                error(['Unknown option (',varargin{opt},').']);
        end
    end
end

%fields{1} = 'rock_group';
fields = {};
if any(strcmpi(fieldsets,'major'))
    fields = {fields{1:end},major{1:end}};
end
if any(strcmpi(fieldsets,'trace'))
    fields = {fields{1:end},trace{1:end}};
end
if any(strcmpi(fieldsets,'ree'))
    fields = {fields{1:end},ree{1:end}};
end
if any(strcmpi(fieldsets,'index'))
    fields = {fields{1:end},index{1:end}};
end

% produces training and testing datasets
metafields = {'rock_group','rock_type','sample_id', ...
    'author','title','journal','year','doi','bibtex'};
nmeta = length(metafields);
traindata = data( :, ...
    {metafields{1:end}, ...
    fields{1:end}} );

% assign serial index values to help identify
traindata.index = [1:height(traindata)]';

% to look at metaigneous and metasedimentary
ind = rockgroup(data,'metasedimentary');
traindata.rock_group(ind) = {'metasedimentary'};
ind = rockgroup(data,'metaigneous');
traindata.rock_group(ind) = {'metaigneous'};

traindata.rock_group = lower(traindata.rock_group);

% treat values that are NaN or < 0
for i = 1:length(fields)
    if usenan
        % remove fields with values < 0
        ind = traindata{:,fields{i}} < 0;
        traindata = traindata(~ind,:);
    else
        % change values that are below detection to NaN
        ind = traindata{:,fields{i}} > 0;
        traindata = traindata(ind,:);
    end
end

% subset data
ind = sum(isnan(traindata{:,nmeta+1:end}),2) == 0;
traindata = traindata(ind,:);

[ind_ig,ind_sed,ind_metaig,ind_metased] = rgindex(traindata.rock_group);

switch option
    case 1
        ind = ind_ig | ind_sed;
    case 2
        ind = ind_metig | ind_metased;
    case 3 % default
        % make a plot of Mahalanobis distance to illustrate metamorphic
        % rocks are similar distributions to their unmetamorphosed
        % counterparts.
        
        [dsquared_ig,mu,V] = mahalanobis(traindata{ind_ig,fields},centering_fcn);
        dsquared_metaig = mahalanobis(traindata{ind_metaig,fields},mu,V,centering_fcn);

        [dsquared_sed,mu,V] = mahalanobis(traindata{ind_sed,fields},centering_fcn);
        dsquared_metased = mahalanobis(traindata{ind_metased,fields},mu,V,centering_fcn);

        figure;
        subplot(221);
        histogram(dsquared_ig,'Normalization','probability', ...
            'BinEdges',[0:2:200],'DisplayStyle','stairs');
        hold on;
        histogram(dsquared_metaig,'Normalization','probability', ...
            'BinEdges',[0:2:200],'DisplayStyle','stairs');
        xlabel('Mahalanobis Distance (D^2) to igneous centroid');
        legend('igneous','metaigneous');
        xlim([0 60]);
        ylabel('PDF');

        subplot(222);
        histogram(dsquared_sed,'Normalization','probability', ...
            'BinEdges',[0:2:200],'DisplayStyle','stairs');
        hold on;
        histogram(dsquared_metased,'Normalization','probability', ...
            'BinEdges',[0:2:200],'DisplayStyle','stairs');
        xlabel('Mahalanobis Distance (D^2) to sedimentary centroid');
        legend('igneous','sedimentary');
        xlim([0 60]);
        ylabel('PDF');

        [dsquared_ig,mu,V] = mahalanobis(traindata{ind_ig | ind_metaig,nmeta+1:end},centering_fcn);
        dsquared_sed = mahalanobis(traindata{ind_sed | ind_metased,nmeta+1:end},mu,V,centering_fcn);

        subplot(212);
        histogram(dsquared_ig,'Normalization','cdf', ...
            'BinEdges',[0:2:200],'DisplayStyle','stairs');
        hold on;
        histogram(dsquared_sed,'Normalization','cdf', ...
            'BinEdges',[0:2:200],'DisplayStyle','stairs');
        xlabel('Mahalanobis Distance (D^2) to igneous centroid');
        legend('meta+igneous','meta+sedimentary');
        xlim([0 200]);
        ylabel('CDF');

        % change metamorphic rock groups to igneous and sedimentary
        traindata.rock_group(ind_metaig) = {'igneous'};
        traindata.rock_group(ind_metased) = {'sedimentary'};

        ind = ind_ig | ind_sed | ind_metaig | ind_metased;
end

% filter for final dataset to use for training and testing
traindata = traindata(ind,:);

% split data into a training and testing dataset
N = height(traindata);
ftrain = (1 - reserve);
if equal
    iig = find(rockgroup(traindata,'igneous'));
    ised = find(rockgroup(traindata,'sedimentary'));

    Ni = length(iig);
    Ns = length(ised);
    
    if Ni > Ns
        N = Ns;
    else
        N = Ni;
    end
    
    iind = sort(randperm(Ni,ceil(ftrain*N)));
    sind = sort(randperm(Ns,ceil(ftrain*N)));
    
    ind = false([Ns + Ni,1]);
    
    ind(sort([iig(iind); ised(sind)])) = 1;    
else
    tind = sort(randperm(N,ceil(ftrain*N)));

    ind = false([N,1]);
    ind(tind) = 1;
end

testdata = traindata(~ind,:);
traindata = traindata(ind,:);

[ind_ig,ind_sed,ind_metaig,ind_metased] = rgindex(traindata.rock_group);

% compute pca (including recentered data if necessary
[~,coeff, scores, latent, tsquared, explained, mu, scale_ig] = db_pca(traindata(ind_ig,:),'class','rock_type',centering_fcn);

[~,coeff, scores, latent, tsquared, explained, mu, scale_sed] = db_pca(traindata(ind_sed,:),'class','rock_type',centering_fcn);


if ~strcmp(centering_fcn,'ilr')
    [traindata{:,fields}, coeff, scores, latent, tsquared, explained, mu] = ...
        db_pca(traindata,'class','rock_type',centering_fcn);
else
    for i = 1:length(fields)-1
        cfields{i} = ['coord_',num2str(i)];
    end
    [traindata{:,cfields}, coeff, scores, latent, tsquared, explained, mu, scale_all] = ...
        db_pca(traindata,'class','rock_type',centering_fcn);
end

if pcaflag
    keep = quantile(tsquared,0.95);

    traindata = traindata(tsquared < keep,:);
end

if strcmp(centering_fcn,'clr')
    testdata{:,fields} = clr(testdata{:,fields});
elseif strcmp(centering_fcn,'ilr')
    testdata{:,cfields} = clr(testdata{:,fields})*scale_all;
end

return


function [ind_ig,ind_sed,ind_metaig,ind_metased] = rgindex(rg)

ind_ig = strcmpi(rg,'igneous');
ind_sed = strcmpi(rg,'sedimentary');

ind_metaig =  strcmpi(rg,'metaigneous');
ind_metased = strcmpi(rg,'metasedimentary');

return