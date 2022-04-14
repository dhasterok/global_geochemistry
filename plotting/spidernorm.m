function [t,p,xref,spiderdata] = spidernorm(spiderdata,varargin)
% SPIDERNORM - computes trace element spider diagrams
% 
%   Options:
%       'Elements'      List of elements. The accompanying list should be
%                       in the order for plotting on a spider diagram.  The
%                       argument needs to be a cell array of element char
%                       codes.k
%
%       'RefData'       Reference table with data for normalizing
%                       one spiecies versus another
%
%       'RefField'      Field used for normalizing geochemistry, default
%                       'sio2' if refdata are included.  Requires an
%                       additional argument with field name.
%
%       'RefVal'        Value used for normalizing geochemistry, requires a
%                       value of RefSpecies in the same units as the data.
%
%       'Logfile'       Write a log file?  Requires a filename as an
%                       additional argument.
%
%       'NormRef'       Reference geochemistry for normalizing spider plot.
%                       Requires an additional string argument with code
%                       for normalizing material.  Default is 'G13', N-MORB
%                       by Gale et al. (2013).
%
%       'Layer'         Some references have multiple materials with
%                       concentrations for normalizing spider plot.
%                       Include layer to indicate which reference material
%                       to use (see earthref.xlsx).
%
%       'PlotSigma'     Adds specified n sigma to compatibility diagram,
%                       use positive values for sigma and negative for
%                       standard error
%
%       'SigmaStyle'    Style of plot for uncertainties, either as a 
                
addpath ref_models

% default values
normref = 'G13';
reffield = '';
xref = [];
sigmaflag = 0;
sigmastyle = 'fill';
ellist = {'Cs', 'Rb', 'Ba', 'Th', ...
    'U', 'Nb', 'Ta', 'K', ...
    'La', 'Ce', 'Pb', 'Mo', ...
    'Pr', 'Sr', 'Ga', 'Zr', 'Hf', ...
    'Nd', 'Sm', 'Eu', 'Li', ...
    'Ti', 'Gd', 'Tb', 'Dy', ...
    'Y', 'Ho', 'Er', ...
    'Tm', 'Yb', 'Lu', 'Zn', ...
    'Mn', 'V', 'Sc', 'Co', ...
    'Cu', 'Ni', 'Cr'};

Q = [0.0001 0.5 2.5 10:10:90 97.5 99.5 99.999]'/100;

opt = 1;
logfile = 0;
if nargin > 1
    while opt <= nargin - 1
        switch lower(varargin{opt})
            case 'refdata'
                refdata = varargin{opt+1};
                if ~isempty(reffield)
                    reffield = 'sio2';
                end
                if ~isempty(xref)
                    xref = 'median';
                end
                opt = opt + 2;
            case 'elements'
                ellist = varargin{opt+1};
                opt = opt + 2;
            case 'reffield'
                reffield = varargin{opt+1};
                opt = opt + 2;
            case 'refval'
                xref = varargin{opt+1};
                opt = opt + 2;
            case 'logfile'
                fid = fopen('spidernorm.log','a+');
                logfile = 1;
                opt = opt + 2;
            case 'normref'
                normref = varargin{opt+1};
                opt = opt + 2;
            case 'layer'
                layer = varargin{opt+1};
                opt = opt + 2;
            case 'plotsigma'
                nsigma = varargin{opt+1};
                opt = opt + 2;
            case 'sigmastyle'
                sigmastyle = varargin{opt+1};
                opt = opt + 2;
            otherwise
                error(['Unknown option, ',varargin{opt}]);
        end
    end
end

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

% make sure all oxide data are in ppm format
[elfield,spiderdata] = ox2ppm(ellist,spiderdata);

if ~isempty(reffield);
    [~,refdata] = ox2ppm(ellist,refdata);

    % reference field for adjusting the geochemistry to global distribution
    xv = refdata{:,reffield};
    x = spiderdata{:,reffield};
    if isempty(xref) | strcmp(xref,'median')
        xref = median(spiderdata{x>0,reffield});
        if logfile
            fprintf(fid,'  %s: %.2f %% median\n',reffield,xref);
        end
    end
end

for i = 1:length(ellist)
    % data for spider
    y = spiderdata{:,elfield{i}};
    
    if ~isempty(reffield) % normalize to reference composition
        % reference data
        yv = refdata{:,elfield{i}};
        
        % accounts for censored data
        indv = ~isnan(xv) & ~isnan(yv);
        ind = ~isnan(x) & ~isnan(y);

        [yref,Qy,ymodel] = adjbyquantile2(xv(indv),yv(indv),xref,x(ind),y(ind));
        

        %if sum(ind) > 0
        %    spiderdata{ind,[elfield{i},'_adj']} = yref;
        %end
        % ignores uncensored data
        %ind = xv > 0 & yv > 0;
        %ind2 = x > 0 & y > 0;
        %[yref,Qy] = adjbyquantile(xv(ind),yv(ind),xref,x(ind2),y(ind2));
        %el(i,1) = nanmedian(log10(yref));
        %el(i,2) = nanstd(log10(yref));

    else % no normalization
        ind = ~isnan(y);
        
        ymodel = gausscensor(y(ind),'scale','log','quantiles',Q);
    end
    
    el(i,1) = ymodel.mu/log(10);
    el(i,2) = ymodel.sigma/log(10);
    el(i,3) = sum(ind);
    
    if logfile
        fprintf(fid,'  %s: %.2f %% censored\n',ellist{i},sum(y <= 0)/(sum(ind) + sum(y <= 0))*100);
    end
        
end

% normalize element data by reference values
el_norm(:,1) = log10(10.^el(:,1)./refchem{1,elfield}');
el_norm(:,2) = el(:,2)./log10(refchem{1,elfield}');

% put results into a table
t.N = el(:,3);
t.mu = el(:,1);
t.sigma = el(:,2);
t.mu_norm = el_norm(:,1);
t.sigma_norm = el_norm(:,2);

t = struct2table(t);
t.Properties.RowNames = ellist';
t.Properties.DimensionNames = {'Element', 'Concentration Data'};
t.Properties.VariableUnits = {'' 'log10 (ppm)' 'log10 (ppm)' 'log10 (ppm)' 'log10 (ppm)'};

% add reference adjustment information to table
t = addprop(t,'Adjusted','table');
t = addprop(t,'RefField','table');
t = addprop(t,'RefValue','table');

if isempty(reffield)
    t.Properties.CustomProperties.Adjusted = 'no';
else
    t.Properties.CustomProperties.Adjusted = 'yes';
end
t.Properties.CustomProperties.RefField = reffield;
t.Properties.CustomProperties.RefValue = xref;

% uncomment this to check and see if you are a moron
%[[1:27]' log10(el(:,1)) log10(refchem{1,elfield}') el_norm]

elid = [1:length(ellist)];

%figure;
if nsigma ~= 0
    ind = ~isnan(t.mu_norm(:,1));
    x = [elid(ind); elid(ind)];
    if nsigma > 0
        y = [t.mu_norm(ind,1) t.mu_norm(ind,1)]' + nsigma*[-t.sigma_norm(ind) t.sigma_norm(ind)]';
    else
        y = [t.mu_norm(ind,1) t.mu_norm(ind,1)]' + nsigma*[-t.sigma_norm(ind)./sqrt(t.N(ind)) t.sigma_norm(ind)./sqrt(t.N(ind))]';
    end

    switch sigmastyle
        case 'errorbar'
            plot(x,y,'Color',[0.5 0.5 0.5]);
            hold on;
            scatter(elid,t.mu_norm(:,1),'o','filled');
        case 'fill'
            ind = isnan(t.mu_norm(:,1));
            x(2,:) = fliplr(x(2,:));
            y(2,:) = fliplr(y(2,:));
            x = x';
            y = y';
            fill(x(:),y(:),[0.8 0.8 0.8],'LineStyle','none');
            hold on;
        otherwise
            error('Unknown SigmaStyle');
    end
end
p = plot(elid,t.mu_norm(:,1),'-');

set(gca,'XTick',elid,'XTickLabel',ellist);
xlim([elid(1)-1 elid(end)+1]);

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