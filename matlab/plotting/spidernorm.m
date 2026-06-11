function [t,h,refval,data] = spidernorm(varargin)
% SPIDERNORM - computes trace element spider diagrams
% 
%   [t,p,xref,spiderdata] = spidernorm(spiderdata,...) computes a trace
%   element compatibility diagram from the table of data spider data.  The
%   function returns a table of normalized values (t) and the relevant
%   statistics, the figure handle (h), the reference value (refval) for
%   normalizing the geochemistry (mean, median, or user specified value),
%   and the data used for normalization.
%
%   Options/value pairs:
%       'Axes'          Axes object for plotting, default is gca
%
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
%       'Quantiles'     Quantile levels (0,1) for plotting, sets 'Style'
%                       automatically to 'Quanta'.  Default quantiles are
%                       [0.05 0.25 0.5 0.75 0.95].
%
%       'Style'         Style 'MeanSD' plots mean and standard
%                       deviation; 'MeanSE' plots mean and standard
%                       error; 'Quanta' (default) plots quantiles.
%
%       'Color'         Color matrix for plotting, each row must be a color
%                       triplet [0,1]
                
addpath ref_models

% default values
sigmaflag = 0;
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

%Q = [0.0001 0.5 2.5 10:10:90 97.5 99.5 99.999]'/100;
Q = [0.05 0.25 0.5 0.75 0.95];

% parse inputs
p = inputParser;

addRequired(p,'DataTable',@istable);
addParameter(p,'Axes',[],@isgraphics);
addParameter(p,'NormRef','G13',@ischar);
addParameter(p,'Layer','',@ischar);
addParameter(p,'RefData',[],@istable);
addParameter(p,'RefField','',@ischar);
addParameter(p,'RefVal',[],@isnumeric);
addParameter(p,'Elements',ellist,@iscellstr);
addParameter(p,'LogFile',false,@islogical);
addParameter(p,'Quantiles',Q,@isnumeric);
addParameter(p,'Color',[0 0.4470 0.7410],@isnumeric);
addParameter(p,'Style','Quanta',@ischar);

parse(p,varargin{:});

data = p.Results.DataTable;
ax = p.Results.Axes;
normref = p.Results.NormRef;
layer = p.Results.Layer;
refdata = p.Results.RefData;
reffield = p.Results.RefField;
refval = p.Results.RefVal;
ellist = p.Results.Elements;
logflag = p.Results.LogFile;
style = p.Results.Style;

Q = p.Results.Quantiles;
C = p.Results.Color;

if isempty(ax)
    ax = gca;
end

% if reference data are included, there needs to be a reference field and
% associated value to normalize the data to.
if ~isempty(refdata)
    if ~isempty(reffield)
        reffield = 'sio2';
    end
    if ~isempty(refval)
        refval = 'median';
    end
end

% start a log file if necessary
if logflag
    fid = fopen('spidernorm.log','a+');
end

% reference for normalizing spider diagram
eref = readtable('earthref.xlsx');

% find indexes of reference data
refind = (strcmpi(eref.model,normref) | strcmpi(eref.reference,normref)) & eref.sigma == 0;

% if multiple layers, filter for desired layer
if sum(refind) > 1
    if isempty(layer)
        error('Add layer for reference.');
    end
    refind = refind & strcmpi(eref.layer,layer);
else
    layer = eref.layer{refind};
end

% reference for normalizing
refchem = eref(refind,:);

% make sure all oxide data are in ppm format
if any(contains(data.Properties.VariableNames,'ppm'))
    [elfield,data] = ox2ppm(ellist,data);
else
    % when the data are in ppm, but there is no '_ppm' appended to field
    % names
    elfield = lower(ellist);
    data.Properties.VariableNames = lower(data.Properties.VariableNames);
    if ~isempty(refdata)
        refdata.Properties.VariableNames = lower(refdata.Properties.VariableNames);
    end
    refchem.Properties.VariableNames = regexprep(refchem.Properties.VariableNames,'_ppm','');
end

if ~isempty(reffield)
    [~,refdata] = ox2ppm(ellist,refdata);

    % reference chemical field for adjusting the geochemistry to global distribution
    xv = refdata{:,reffield};
    x = data{:,reffield};
    if isempty(refval) || strcmp(refval,'median')
        refval = median(data{x>0,reffield});
        if logflag
            fprintf(fid,'  %s: %.2f %% median\n',reffield,refval);
        end
    end
end

for i = 1:length(ellist)
    % data for spider
    y = data{:,elfield{i}};
    
    if ~isempty(reffield) % normalize to reference composition
        %if i == 1
        %    disp(['Normalization to ',elfield{i},'...']);
        %end
        % reference data
        yv = refdata{:,elfield{i}};
        
        % accounts for censored data
        indv = ~isnan(xv) & ~isnan(yv);
        ind = ~isnan(x) & ~isnan(y);

        [yref,yq,ymodel] = adjbyquantile2(xv(indv),yv(indv),refval,x(ind),y(ind));

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
        %if i == 1
        %    disp('No normalization to species...');
        %    figure
        %    histogram(log10(data.cs_ppm));
        %    hold on;
        %end
        ind = ~isnan(y);

        % censored data?
        if all(y > 0)       % no
            %disp([ellist{i},' uncensored'])
            switch style
                case {'MeanSD','MeanSE'}
                    ymodel.mu = mean(log10(y(ind)));
                    ymodel.sigma = std(log10(y(ind)));
                case 'Quanta'
                    yq = quantile(log10(y),Q);
                    %if i == 1
                    %    yl = get(gca,'YLim');
                    %    plot(yq(2)*[1 1],yl);
                    %end
            end
        else                % yes
            %disp([ellist{i},' censored'])
            [ymodel,yq] = gausscensor(y(ind),'Scale','log10','Quantiles',Q);
            %if i == 1
            %    yl = get(gca,'YLim');
            %    plot(yq(2,2)*[1 1],yl);
            %end
        end
    end
    
    switch style
        case {'MeanSD','MeanSE'}
            el(i,1) = ymodel.mu/log(10);
            el(i,2) = ymodel.sigma/log(10);
            el(i,3) = sum(ind);
        case 'Quanta'
            if width(yq) == 2
                el(i,:) = [yq(:,2)' sum(ind)];
            else
                el(i,:) = [yq sum(ind)];
            end
    end
    
    if logflag
        fprintf(fid,'  %s: %.2f %% censored\n',ellist{i},sum(y <= 0)/(sum(ind) + sum(y <= 0))*100);
    end
        
end

% put results into a table
switch style
    case {'MeanSD','MeanSE'}
        % normalize element data by reference values
        el_norm(:,1) = log10(10.^el(:,1)./refchem{1,elfield}');
        el_norm(:,2) = el(:,2)./log10(refchem{1,elfield}');

        t.N = el(:,end);
        t.mu = el(:,1);
        t.mu_norm = el_norm(:,1);

        switch style
            case {'MeanSD'}
                t.sigma = el(:,2);
                t.sigma_norm = el_norm(:,2);
            case 'MeanSE'
                t.se = el(:,2)./sqrt(t.N);
                t.se_norm = el_norm(:,2);
        end

        t = struct2table(t);
        t.Properties.RowNames = ellist';
        t.Properties.DimensionNames = {'Element', 'Concentration Data'};
        t.Properties.VariableUnits = {'' 'log10 (ppm)' 'log10 (ppm)' 'log10 (ppm)' 'log10 (ppm)'};
        
    case 'Quanta'
        el_norm = log10(10.^el(:,1:end-1)./repmat(refchem{1,elfield}',1,width(el)-1));
        for i = 1:length(Q)
            qtxt{i} = ['Q',num2str(Q(i)*100)];
        end
        qtxt{end+1} = 'Ndata';
        %t = array2table([el(:,end) el(:,1:end-1) el_norm],'VariableNames',qtxt,'RowNames',ellist');
        t = array2table([el_norm el(:,end)],'VariableNames',qtxt,'RowNames',ellist');
end

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
t.Properties.CustomProperties.RefValue = refval;

% uncomment this to check and see if you are a moron
%[[1:27]' log10(el(:,1)) log10(refchem{1,elfield}') el_norm]

elid = [1:length(ellist)];

%figure;
x = [elid; elid];
switch style
    case 'MeanSD'
        ind = ~isnan(t.mu_norm(:,1));
        y1 = [t.mu_norm(ind,1) t.mu_norm(ind,1)]' + [-t.sigma_norm(ind) t.sigma_norm(ind)]';
        y2 = [t.mu_norm(ind,1) t.mu_norm(ind,1)]' + 2*[-t.sigma_norm(ind) t.sigma_norm(ind)]';

        h(3) = fillinterval(ax,x,y2,C,0.15); hold on;
        h(2) = fillinterval(ax,x,y1,C,0.3);
        h(1) = plot(ax,elid,t.mu_norm(ind,1),'Color',C,'LineWidth',0.75)
        yl = [floor(min(y2(1,:),[],'omitnan')) ceil(max(y2(2,:),[],'omitnan'))];
    case 'MeanSE'
        ind = ~isnan(t.mu_norm(:,1));
        y1 = [t.mu_norm(ind,1) t.mu_norm(ind,1)]' + [-t.se_norm(ind) t.se_norm(ind)]';
        y2 = [t.mu_norm(ind,1) t.mu_norm(ind,1)]' + 2*[-t.se_norm(ind) t.se_norm(ind)]';

        h(3) = fillinterval(ax,x,y2,C,0.15);
        h(2) = fillinterval(ax,x,y1,C,0.3);
        h(1) = plot(ax,elid,t.mu_norm(ind,1),'Color',C,'LineWidth',0.75)
        yl = [floor(min(y2(1,:),[],'omitnan')) ceil(max(y2(2,:),[],'omitnan'))];
    case 'Quanta'
        switch length(Q)
            case 1
                ind = ~isnan(t{:,1});
                h = plot(ax,elid(ind),t{ind,1}','Color',C,'LineWidth',0.75);
                yl = [floor(min(t{:,1},[],'omitnan')) ceil(max(t{:,1},[],'omitnan'))];
            case 2
                fillinterval(ax,x,t{:,[1 2]}',C,0.3);
                yl = [floor(min(t{:,1},[],'omitnan')) ceil(max(t{:,2},[],'omitnan'))];
            case 3
                ind = ~isnan(t{:,2});
                h(2) = fillinterval(ax,x,t{:,[1 3]}',C,0.3);
                h(1) = plot(ax,elid(ind),t{ind,2}','Color',C,'LineWidth',0.75);
                yl = [floor(min(t{:,1},[],'omitnan')) ceil(max(t{:,3},[],'omitnan'))];
            case 5
                ind = ~isnan(t{:,3});
                h(3) = plot(ax,x(ind),[t{ind,[1 5]}]','Color',C,'LineWidth',0.25);
                h(2) = fillinterval(ax,x,t{:,[2 4]}',C,0.3);
                h(1) = plot(ax,elid(ind),t{ind,3}','Color',C,'LineWidth',0.75);
                yl = [floor(min(t{:,1},[],'omitnan')) ceil(max(t{:,5},[],'omitnan'))];
            otherwise
                for i = 1:length(Q)
                    plot(ax,elid,t{:,i},'Color',C);
                end
        end
end

% switch sigmastyle
%     case 'errorbar'
%         plot(ax,x,y,'Color',[0.5 0.5 0.5]);
%         hold on;
%         scatter(ax,elid,t.mu_norm(:,1),'o','filled');
%     case 'fill'
%         ind = isnan(t.mu_norm(:,1));
%         x(2,:) = fliplr(x(2,:));
%         y(2,:) = fliplr(y(2,:));
%         x = x';
%         y = y';
%         fill(ax,x(:),y(:),[0.8 0.8 0.8],'LineStyle','none');
%         hold on;
%     otherwise
%         error('Unknown SigmaStyle');
% end
%p = plot(ax,elid,t.mu_norm(:,1),'-');

%set(ax,'XTick',elid,'XTickLabel',ellist);
ax.XTick = elid;
ax.XTickLabel = ellist;
ax.XLim = [elid(1)-1 elid(end)+1];
logax(yl,'Axes',ax,'Axis','y')

if logflag
    fclose(fid);
end

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


function h = fillinterval(ax,x,y,C,alpha)

%ind = isnan(y(1,:));
x(2,:) = fliplr(x(2,:));
y(2,:) = fliplr(y(2,:));

x = x';
y = y';

h = fill(ax,x(:),y(:),C,'LineStyle','none','FaceAlpha',alpha);

return