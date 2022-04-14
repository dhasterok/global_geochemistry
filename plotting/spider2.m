function t = spider(data,varargin)
% SPIDER - computes 
%
%   Options:
%       'Elements'      list of elements in cell array
%
%       'NormRef'       normalizing reference material found in
%                       earthref.xlsx including the layer separated by a
%                       colon (model:layer), e.g., 'RG03:layer', the
%                       default is NMORB as defined by G13, Gale et al.
%                       (2013)
%
%                       use 'mean' for the mean of all subset values, and
%                       'global' for the median of all the data
%
%       'TraceNorm'     normalize trace elements to a chemical species,
%                       include in a cell array both the species and the
%                       value to normalize to e.g., {'sio2', 62}
%
%                       if TraceNorm is not provided or value is 'none', the
%                       data will not be normalized
%
%       'Subset'        subset indices; to include multiple subsets use a
%                       logical array or specify 'field' and provide the
%                       additional option of 'BinEdges' for a numeric
%                       field or 'SubsetField' for a text field
%
%       'SubsetField'   field for subsetting data into bins or by unique
%                       values
%
%       'BinEdges'      numerical field with edges of bins for producing
%                       subsets
%
%       'Subtext'       text for subsets to be put in legend, can also be
%                       used choose selected subsets for plotting
%
%       'Color'         color for spider plot, for multiple subsets
%                       specify a color for each
%
%       'LogFile'       Produces a log in spidernorm.log

addpath ref_models

% default values
normref = 'G13';
reffield = '';
xref = 'median';
ellist = {'Cs', 'Rb', 'Ba', 'Th', ...
    'U', 'Nb', 'Ta', 'P', 'K', ...
    'La', 'Ce', 'Pb', 'Mo', ...
    'Pr', 'Sr', 'Ga', 'Zr', 'Hf', ...
    'Nd', 'Sm', 'Eu', 'Li', ...
    'Ti', 'Gd', 'Tb', 'Dy', ...
    'Y', 'Ho', 'Er', ...
    'Tm', 'Yb', 'Lu', 'Zn', ...
    'Mn', 'V', 'Sc', 'Co', ...
    'Cu', 'Ni', 'Cr'};
subset = 'all';
subtext = '';
colour = [];
edges = [];

opt = 1;
logfile = 0;
if nargin >= opt + 1
    while opt + 1 < nargin
        switch lower(varargin{opt})
            case 'elements'
                ellist = varargin{opt+1};
                opt = opt + 2;
            case 'tracenorm'
                reffield = varargin{opt+1}{1};
                xref = varargin{opt+1}{2};
                opt = opt + 2;
            case 'normref'
                [normref,layer] = strtok(varargin{opt+1},':');
                if ~isempty(layer)
                    layer = layer(2:end);
                end
                opt = opt + 2;
            case 'layer'
                layer = varargin{opt+1};
                opt = opt + 2;
            case 'subset'
                subset = varargin{opt+1};
                opt = opt + 2;
            case 'subsetfield'
                field = varargin{opt+1};
                opt = opt + 2;
            case 'binedges'
                edges = varargin{opt+1};
                opt = opt + 2;
            case 'subtext'
                subtext = varargin{opt+1};
                opt = opt + 2;
            case 'color'
                colour = varargin{opt+1};
                opt = opt + 2;
            case 'logfile'
                fid = fopen('spidernorm.log','a+');
                logfile = 1;
                opt = opt + 1;
            otherwise
                error(['Unknown option, ',varargin{opt}]);
        end
    end
end

if islogical(subset)
    % do nothing, data are already subset
elseif strcmp(subset,'field')
    if isempty(edges)
        % if subset chars/strings are not provided
        if isempty(subtext)
            subtext = unique(data{:,field})
        end
        % indices for char/string bins
        subset = logical(zeros([height(data),length(subtext)]));
        for i = 1:length(subtext)
            subset(:,i) = strcmp(data{:,field},subtext{i});
        end
    else
        % indicies for numeric bins
        subset = logical(zeros([height(data),length(edges)-1]));
        for i = 1:length(edges)-1
            subset(:,i) = edges(i) <= data{:,field} & data{:,field} < edges(i+1);
            subtext{i} = [num2str(edges(i)),' to ',num2str(edges(i+1))];
        end
        if isempty(colour)
            colour = flipud([flipud(autumn(ceil((length(edges)-1)/2))); ...
                winter(ceil((length(edges)-1)/2))]);
        end
    end 
elseif strcmp(subset,'all')
    subset = true(height(data),1);
    subtext = {'all'};
end

if ~strcmp(normref,{'mean','global'})
    % reference for normalizing spider diagram
    eref = readtable('earthref.xlsx');

    % find indexes of reference data
    ind = (strcmpi(eref.model,normref) | strcmpi(eref.reference,normref)) & eref.sigma == 0;

    % if multiple layers, filter for desired layer
    if sum(ind) > 1
        if isempty(layer)
            error('Add layer for reference.');
        end
        ind = ind & strcmpi(eref.layer,layer);
    else
        layer = eref.layer{ind};
    end

    % reference for normalizing
    refchem = eref(ind,:);
end

% make sure all oxide data are in ppm format
[elfield,data] = ox2ppm(ellist,data);

% normalize chemistry if necessary
if ~isempty(reffield) & ~strcmp(reffield,'none')
    if isempty(xref) | strcmp(xref,'median')
        xref = median(data{:,reffield});
    end
    if logfile
        fprintf(fid,'  %s: %.2f %% median\n',reffield,xref);
    end
    [smooth_model,adjdata] = tracenorm(data,'refspecies',reffield,'refval',xref,'elements',ellist);
end

if ~iscolumn(subtext)
    subtext = subtext';
end
    
t = table(subtext,'VariableNames',{'subset'});
for i = 1:size(subset,2)        % each subset
    for j = 1:length(ellist)    % each element
        if ~isempty(reffield)
            y = adjdata{subset(:,i),elfield{j}};
        else
            y = data{subset(:,i),elfield{j}};
        end
        
        % unnormalized by reference chemical species
        ymodel = gausscensor(y,'scale','log');
        
%         if j == 1
%             histogram(log(adjdata{adjdata{:,elfield{j}}>0,elfield{j}}),'Normalization','probability')
%             hold on;
%             histogram(log(y(y>0)),'Normalization','probability');
%             
%             plot(ymodel.mu*[1 1],get(gca,'YLim'),'k-');
%             plot(log(refchem{1,elfield{1}})*[1 1],get(gca,'YLim'),'r-');
%             title([elfield{j},' ',subtext{i}]);
%             pause
%             hold off;
%         end
        
        elfieldm{j} = [elfield{j},'_mu'];
        elfieldmn{j} = [elfield{j},'_mu_norm'];
        elfields{j} = [elfield{j},'_sigma'];
        elfieldsn{j} = [elfield{j},'_sigma_norm'];
        warning off;
        t{i,elfieldm{j}} = ymodel.mu;
        t{i,elfields{j}} = ymodel.sigma;

        ind = ~isnan(y);
        t{i,[elfield{j},'_N']} = sum(ind);
        
        if logfile
            fprintf(fid,'  %s: %.2f %% censored\n',ellist{j},sum(y <= 0)/(sum(ind) + sum(y <= 0))*100);
        end
    end

    % uncomment this to check and see if you are a moron
    %[[1:27]' log10(el(:,1)) log10(refchem{1,elfield}') el_norm]
end

% normalize element data by reference values
switch normref
    case 'mean'
        ref = repmat(nanmean(exp(t{:,elfieldm})),size(subset,2),1);
    case 'global'
        ref = repmat(nanmedian(adjdata{:,elfield}),size(subset,2),1);
    otherwise
        ref = repmat(refchem{1,elfield},size(subset,2),1);
end

t{:,elfieldmn} = log10(exp(t{:,elfieldm})./ref);
t{:,elfieldsn} = (t{:,elfields}/log(10))./log10(ref);

elid = [1:length(ellist)];
if ~isempty(colour)
    for i = 1:size(subset,2)
        p = plot(elid',t{i,elfieldmn}','-');
        hold on;
        set(p,'Color',colour(i,:));
    end
else
    p = plot(elid',t{:,elfieldmn}','-');
end
set(gca,'Box','on','XTick',elid,'XTickLabel',ellist);
xlim([elid(1)-1 elid(end)+1]);
yl = get(gca,'YLim');
hpax([floor(yl(1)) ceil(yl(2))],'y');
switch normref
    case 'mean'
        ylabel(['Abundance/Mean']);
    otherwise
        ylabel(['Abundance/',refchem.layer{1},' (',refchem.model{1},')']);
end
golden        

if ~isempty(subtext)
    legend(subtext,'Location','eastoutside');
end

if logfile
    fclose(fid);
end
%pause;
%close all;

return


% convert oxides to ppm
function [elfield,data] = ox2ppm(ellist,data)

for i = 1:length(ellist)
    elfield{i} = [lower(ellist{i}),'_ppm'];
    
    switch elfield{i}
        case 'k_ppm'
            if any(strcmp(data.Properties.VariableNames,'k2o'))
                data.k_ppm = 2*molecularwt('K')/molecularwt('K2O')*data.k2o*10000;
            end
        case 'na_ppm'
            if any(strcmp(data.Properties.VariableNames,'na2o'))
                data.na_ppm = 2*molecularwt('Na')/molecularwt('Na2O')*data.na2o*10000;
            end    
        case 'ti_ppm'
            if any(strcmp(data.Properties.VariableNames,'tio2'))
                data.ti_ppm = molecularwt('Ti')/molecularwt('TiO2')*data.tio2*10000;
            end
        case 'al_ppm'
            if any(strcmp(data.Properties.VariableNames,'al2o3'))
                data.al_ppm = 2*molecularwt('Al')/molecularwt('Al2O3')*data.al2o3*10000;
            end
        case 'ba_ppm'
            if any(strcmp(data.Properties.VariableNames,'bao'))
                data.ba_ppm = molecularwt('Ba')/molecularwt('BaO')*data.bao*10000;
            end
        case 'fe_ppm'
            if any(strcmp(data.Properties.VariableNames,'feo_tot'))
                data.fe_ppm = molecularwt('Fe')/molecularwt('FeO')*data.feo_tot*10000;
            end    
        case 'ni_ppm'
            if any(strcmp(data.Properties.VariableNames,'nio'))
                data.ni_ppm = molecularwt('Ni')/molecularwt('NiO')*data.nio*10000;
            end
        case 'mn_ppm'
            if any(strcmp(data.Properties.VariableNames,'mno'))
                data.mn_ppm = molecularwt('Mn')/molecularwt('MnO')*data.mno*10000;
            end
        case 'cr_ppm'
            if any(strcmp(data.Properties.VariableNames,'cr2o3'))
                data.cr_ppm = 2/3*molecularwt('Cr')/molecularwt('Cr2O3')*data.cr2o3*10000;
            end
        case 'p_ppm'
            if any(strcmp(data.Properties.VariableNames,'p2o5'))
                data.p_ppm = 2*molecularwt('P')/molecularwt('P2O5')*data.p2o5*10000;
            end
    end
end

return