 function varargout = tracenorm(data,varargin)
% TRACENORM - Normalizes trace element data by quantile.
%
%   smooth_model = tracenorm(data) will compute the scale parameters
%   necessary for producing normalizing trace element data against a
%   reference oxide (default SiO2).
%
%   [smooth_model,adjdata] = tracenorm(data,'refval',value) will normalize
%   the trace element data to the reference species with a concentration of
%   value.  The adjdata will return a table with the same metadata (all
%   non-numeric fields) as data with the normalized trace element
%   concentrations.
%
%   [smooth_model,adjdata] = tracenorm(data,'subset',index,'refval',value)
%   will return the table of normalized data (adjdata) only for the samples
%   data(index,:).
%
%   Not all trace elements that are normalized by default.  Only
%   Cs, Rb, Ba, Th, U, Nb, Ta, K, La, Ce, Pb, Mo, Pr, Sr, Ga, Zr, Hf, Nd,
%   Sm, Eu, Li, Ti, Gd, Tb, Dy, Y, Ho, Er, Tm, Yb, Lu, Zn, Mn, V, Sc, Co, 
%   Cu, Ni, and Cr, are normalized in by default.
%
%   Input additional options as tracenorm(data, 'option', value, ...).
%   User options:
%       'refspecies'    Reference chemical species for determining quantile
%                       levels and normalizing
%       'refval'        Value of refspecies to normalize trace element data
%       'elements'      Change the element list for normalization, input as
%                       a cell array
%       'subset'        Normalize only a subset of the data by providing an
%                       logical array of height(data), where true indicates
%                       the data to be normalized
%       'plot'          Creates a plot or specified plots as given by cell
%                       array (see plot types below)
%
%   Multiple plot types may be specified in a cell array.
%   Plot types:
%       'quantile'      Produces a contour plot of observed and smoothed
%                       quantiles for each trace element as a function of
%                       the reference species 
%       'data'          Plots 2D histograms of the input data for each
%                       trace element
%       'scale'         Produces a plot of log-normal scale factors (mu and
%                       sigma, and a linear least squares fit to the scale
%                       parameters as a function of the reference species
%       'cdfplot'       produces cdf
%       'elquantile'    A histogram of the smoothed quantile levels
%                       associated with each data point
%       'elhist'        Histograms of the trace element concentrations for
%                       the data observed and normalized to the reference
%                       species value
%
%   Note: plots rescale trace element concentrations on a natural log
%   scale.

% defaults
reffield = 'sio2';
xref = [];
ellist = {'Cs', 'Rb', 'Ba', 'Th', ...
    'U', 'Nb', 'Ta', 'K', ...
    'La', 'Ce', 'Pb', 'Mo', ...
    'Pr', 'Sr', 'Ga', 'Zr', 'Hf', ...
    'Nd', 'Sm', 'Eu', 'Li', ...
    'Ti', 'Gd', 'Tb', 'Dy', ...
    'Y', 'Ho', 'Er', ...
    'Tm', 'Yb', 'Lu', 'Zn', ...
    'Mn', 'V', 'Sc', 'Co' ...
    'Cu', 'Ni', 'Cr'};
ind = [];
plot_type = [];

% user-defined options
opt = 1;
if nargin >= opt + 2
    while opt + 1 < nargin
        switch lower(varargin{opt})
            case 'elements'
                ellist = varargin{opt+1};
                opt = opt + 2;
            case 'refspecies'
                reffield = varargin{opt+1};
                opt = opt + 2;
            case 'refval'
                xref = varargin{opt+1};
                opt = opt + 2;
            case 'subset'
                ind = varargin{opt+1};
                opt = opt + 2;
            case 'plot'
                plot_type = varargin{opt+1};
                opt = opt + 2;
            otherwise
                error(['Unknown option, ',varargin{opt}]);
        end
    end
end

% convert oxides from the element list to concentrations in ppm
data = ox2ppm(ellist,data);

% Output results for each element
for i = 1:length(ellist)
    % change element value to same as trace element field
    elfield = [lower(ellist{i}),'_ppm'];
    
    if isempty(ind) % when no subset is provided
        % if no reference is provided, compute the smooth models, no data
        % are adjusted
        if isempty(xref)
            [xq,yq,Yq,model,smooth_model(i)] = adjbyquantile2(data{:,reffield},data{:,elfield});
        else % adjust data to xref
            % initialize adjdata table
            if i == 1
                adjdata = copymetadata(data);
                adjdata{:,reffield} = xref;
            end
            
            % compute normalized data
            [y_at_ref,Qy,ymodel(i),yref,Yref,Q,xq,yq,Yq,model,smooth_model(i)] = adjbyquantile2(data{:,reffield},data{:,elfield},xref);
            %sum(~isnan(data{:,elfield}))
            adjdata{:,elfield} = y_at_ref;
        end
    else % for subset data
        % initialize adjdata table
        if i == 1
            adjdata = copymetadata(data(ind,:));
            adjdata{:,reffield} = xref;
        end
        
        % compute normalized data
        [y_at_ref,Qy,ymodel(i),yref,Yref,Q,xq,yq,Yq,model,smooth_model(i)] = adjbyquantile2(data{:,reffield},data{:,elfield},xref, ...
            data{ind,reffield},data{ind,elfield});
        adjdata{:,elfield} = y_at_ref;
    end

    % produce plots if requested
    if ~isempty(plot_type)
        if i == 1
            Nr(1) = 8;
            Nr(2) = ceil((length(ellist)+1)/Nr(1));
        end
        for j = 1:length(plot_type)
            switch lower(plot_type{j})
                case 'quantile'
                    quantile_plot(Nr,xq,yq,Yq,reffield,ellist{i});
                case 'data'
                    data_plot(Nr,data{:,reffield},data{:,elfield},reffield,ellist{i})
                case 'scale'
                    scale_plot(Nr,model,smooth_model(i),reffield,ellist{i})
                case 'residual'
                    residual_plot(Nr,model,smooth_model(i),reffield,ellist{i})
                case 'cdfplot'
                    cdf_plot(Nr,xq,yq,xref,yref,Yref,Q,reffield,ellist{i})
                case 'elquantile'
                    elquantile(Nr,Qy,ellist{i})
                case 'elhist'
                    elhist(Nr,data{ind,elfield},y_at_ref,ellist{i})
                otherwise
                    error(['Unknown plot type, ',plot_type{j}]);
            end
        end
    end
end

%turn smooth_model structure into a table
smooth_model = struct2table(smooth_model);
smooth_model = addvars(smooth_model,ellist','NewVariableNames','element','Before',1);

% return variables
if nargout > 0
    varargout{1} = smooth_model;
    if ~isempty(xref)
        varargout{2} = adjdata;
    end
end

return


% convert from oxide wt.% to ppm
function data = ox2ppm(ellist,data)

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


% extract only columns with metadata (i.e. non-numeric columns)
function t = copymetadata(t)

ind = false(1,width(t));
for i = 1:width(t)
    if ~strcmp(class(t{:,i}),'double')
        ind(i) = 1;
    end
end

t = t(:,ind);

return


% quantile plots
function quantile_plot(Nr,xq,yq,Yq,reffield,el)

persistent qfig k
try
    figure(qfig);
    k = k + 1;
catch
    qfig = figure;
    k = 1;
end

yl = [floor(min(Yq(:)))-1 ceil(max(Yq(:)))+1];

subplot(Nr(2),Nr(1),k);
hold on;

% 2-d histogram/quantile plots
%plot(xq',yq','r-');
plot(xq',Yq','b-');

title(el);
xlabel([reffield,' (wt.%)']);
ylabel(['log(',el,', ppm)']);

golden;
xlim([42 82]);
ylim([-10 14]);

set(gca,'XTick',[40:5:80],'YTick',[-10:2:14]);
set(gca,'Box','on','TickDir','out');
    
return


% data histograms
function data_plot(Nr,xv,yv,reffield,el)

persistent qfig k
try
    figure(qfig);
    k = k + 1;
catch
    qfig = figure;
    k = 1;
end

Q = quantile(log(yv(yv > 0)),[0.01 99.9]/100);

yl = [floor(Q(1))-1 ceil(Q(2))+1];
Y = linspace(yl(1),yl(2),30);
X = linspace(min(xv),max(xv),21);
n = hist2d(xv(yv > 0),log(yv(yv > 0)),X,Y);
    
subplot(Nr(2),Nr(1),k);
hold on;

imagesc(X,Y,log10(n));
axis xy;

%colormap(flipud(gray));
colormap([1 1 1;parula])
caxis([-0.1 3]);
colorbar;

title(el);
xlabel([reffield,' (wt.%)']);
ylabel(['log(',el,', ppm)']);

golden;
xlim([42 82]);
ylim(yl);
set(gca,'XTick',[40:5:80],'YTick',[-12:2:12]);
set(gca,'Box','on','TickDir','out');
    
return


% mu-sigma plots
function scale_plot(Nr,model,smooth_model,reffield,el)

persistent qfig k
try
    figure(qfig);
    k = k + 1;
    
    % if creating the figure, add a legend
    subplot(Nr(2),Nr(1),Nr(2)*Nr(1));
    s = 4*(log([10 30 100 300 1000 3000])-1);
    scatter(ones([1,length(s)]),[1:length(s)],s,'o','filled');
    golden;
    for i = 1:length(s)
        text(1.1,i,num2str(exp(s(i)/4+1)));
    end
    set(gca,'Visible','off');
catch
    qfig = figure;
    k = 1;
end
subplot(Nr(2),Nr(1),k);
hold on;

dx = model.x(2) - model.x(1);
minx = model.x(1) - dx/2;
maxx = model.x(end) - dx/2;

% mu
t = linspace(minx,maxx,20)';
yyaxis left
%scatter(model.x,mu,18,'o','filled');
%scatter(model.x,model.mu,4*(log(model.N)-1),'o','filled');
plot(t,[ones(size(t)) t t.^2]*smooth_model.mm,'-');

ystr = ['Predicted \mu_{',el,'}'];
ylabel(ystr);
title(el);
xlabel([reffield,' (wt.%)']);
set(gca,'Box','on');
ylim([floor(min(model.mu)) ceil(max(model.mu))]);
golden

% sigma
yyaxis right
r = [minx maxx]';
%scatter(model.x,sigma,18,'o','filled');
%scatter(model.x,model.sigma,4*(log(model.N)-1),'o','filled');
plot(r,[ones(size(r)) r]*smooth_model.ms,'-');

ystr = ['Predicted \sigma_{',el,'}'];
ylabel(ystr);
set(gca,'Box','on');
ylim([0 2.5]);
golden
set(gca,'XTick',[50:10:80],'TickDir','out');
    
return


% residual plots
function residual_plot(Nr,model,smooth_model,reffield,el)

persistent qfig k
try
    figure(qfig);
    k = k + 1;
catch
    qfig = figure;
    k = 1;
end
subplot(Nr(2),Nr(1),k);
hold on;

ind = 48 < model.x & model.x < 78;
% mu
histogram(model.mu(ind) - [ones(size(model.x(ind))) model.x(ind) model.x(ind).^2]*smooth_model.mm,'BinEdges',[-1:0.1:1]);
ystr = ['N bins'];
ylabel(ystr);
xlabel([el,' (ppm)']);

% sigma
histogram(model.sigma(ind) - [ones(size(model.x(ind))) model.x(ind)]*smooth_model.ms,'BinEdges',[-1:0.1:1]);
set(gca,'Box','on');
golden
set(gca,'TickDir','out');
    
return


function cdf_plot(Nr,xq,yq,xref,yref,Yref,Q,reffield,el)

persistent qfig k
try
    figure(qfig);
    k = k + 1;
        
    subplot(Nr(2),Nr(1),Nr(2)*Nr(1));
    hold on;
    plot(0,0,'ro-');
    plot(0,0,'bo-');
    plot(0,0,'go-');
    golden;
    legend('empirical','ln N fit','smoothed');
    set(gca,'Visible','off');
catch
    qfig = figure;
    k = 1;
end

subplot(Nr(2),Nr(1),k);
hold on;

[~,p] = min(abs(xref - xq(1,:)));

plot(yq(:,p),Q,'ro-');
plot(yref,Q,'bo-');
plot(Yref,Q,'go-');

title(['Quantiles at ',num2str(xref),' wt.% ',reffield]);
ylabel('CDF');
xlabel(['log(',el,', ppm)']);
golden
set(gca,'Box','on');

return


function elquantile(Nr,Qy,el)

persistent qfig k
try
    figure(qfig);
    k = k + 1;
catch
    qfig = figure;
    k = 1;
end

subplot(Nr(2),Nr(1),k);
hold on;
histogram(Qy,'BinWidth',0.05);

title(el);
xlabel('Quantiles of data subset');
xlim([0 1]);
set(gca,'Box','on');

return


% figures that illustrate the normalization process
function elhist(Nr,y,y_at_ref,el)

persistent qfig k
try
    figure(qfig);
    k = k + 1;
    
    subplot(Nr(2),Nr(1),Nr(2)*Nr(1));
    hold on;
    plot(0,0,'.');
    plot(0,0,'.');
    golden;
    legend('observed','adjusted');
    set(gca,'Visible','off');
catch
    qfig = figure;
    k = 1;
end

subplot(Nr(2),Nr(1),k);
hold on;

histogram(log(y(y > 0)),'BinWidth',0.25);
histogram(log(y_at_ref(y_at_ref > 0)),'BinWidth',0.25);

set(gca,'Box','on');
xlabel(['log(',el,', ppm)']);

return