function [x,y,n] = gaussmavg_plot(data,el1,el2,dx,nx,varargin)
% GAUSSMAVG_PLOT - produces a moving average plot
%
%   [x,y,n] = gaussmavg_plot(data,el1,el2,dx,nx) produces a moving average
%   for el2 binned by el1, where el1 and el2 are either fields of data, or
%   ratios of fields in data.  To supply a ratio, provide the fields in a
%   2-element cell array with the first entry being the numerator.  The
%   averaging width is dx*(nx + 1) wide, so nx must be even.  If not, nx
%   will automatically be increased by 1 to force the value to be even.
%   The step el for in each windowed average is dx.
%
%   The return value x is location of window centers for el1, y is the
%   quantile values for el2 (0.025 0.25 0.5 0.75 0.975), and n is the
%   number of data in each window.
%
%   Many of the input fields can be simplified, e.g. 'rb_ppm' can be input
%   as 'Rb', and many property fields can be input with the field name as
%   'Density' or 'Heat Production' instead of 'density_model' or
%   'heat_production'.
%
%   gaussmavg_plot(data,el1,el2,dx,nx,'Hist') will add a histogram subplot
%   beneath the moving average plot.
%
%   There are different averaging schemes for the data:
%       gaussmavg_plot(data,el1,el2,dx,nx,'Method',method_name), which
%       include 'cdf', 'ideal', and 'dhuime'.  'cdf' will incorporate
%       censored data (to non-ratio quantities) and estimate quantiles from
%       the empirical distribution.  'ideal' will also incorporate censored
%       data to estimate quantiles from a gaussian model of the empirical
%       distribution.  'nocensor' will produce fit a gaussian to the data
%       without accounting for censored data.
%
%   Most options are supplied as option, value pairs.
%
%   Options:
%       'logx', 'logy'  linear or log10 for the data
%
%       'Method'        'cdf', 'ideal', or 'nocensor' (see description
%                       above)
%
%       'ErrorType'     'q' will plot quantiles (default) or 'qerr'
%                       to plot the quantile errors
%
%       'Color'         add color triplet for plotting
%
%       'Alpha'         add value from 0 to 1 indicating the level of
%                       transparency for the uncertainty bounds
%
%       'Hist'          will add a histogram of the number of data in each
%                       step width (dt)


xscale = 'linear';
yscale = 'linear';
edges = [];
colour = [];
method = 'cdf';
alpha = [];
plothist = 0;
errtype = 'q';
if nargin > 4
    opt = 1;
    while opt + 5 <= nargin
        switch lower(varargin{opt})
            case 'logx'
                xscale = 'log';
                opt = opt + 1;
            case 'logy'
                yscale = 'log';
                opt = opt + 1;
            case 'color'
                colour = varargin{opt+1};
                opt = opt + 2;
            case 'method'
                method = varargin{opt+1};
                opt = opt + 2;
            case 'alpha'
                alpha = varargin{opt+1};
                opt = opt + 2;
            case 'hist'
                plothist = 1;
                opt = opt + 1;
            case 'errortype'
                errtype = varargin{opt+1};
                opt = opt + 2;
            otherwise
                error(['Unknown option (',varargin{opt},').']);
        end
    end
end

% create axis labels before fixing field names
xlbl = getlabel(el1);
ylbl = getlabel(el2);

% fix field names so that they match those in data
el1 = fixfieldnames(data.Properties.VariableNames,el1);
el2 = fixfieldnames(data.Properties.VariableNames,el2);

ind = data{:,el1} > 0;
data = data(ind,:);
if strcmp(xscale,'log')
    x = [floor(min(log10(data{:,el1}))/dx)*dx : dx : ceil(max(log10(data{:,el1}))/dx)*dx]';
else
    x = [(floor(min(data{:,el1})/dx)-nx/2)*dx : dx : (ceil(max(data{:,el1})/dx)+nx/2)*dx]';
    
end

if iscell(el2) & length(el2) == 2
    % to handle censoring it would require dealing with BDL values in the
    % denominator resulting in right censoring of ratios, or what to do
    % with BDL/BDL values?  Best to skip handling them for now.
    ind = data{:,el2{1}} > 0 & data{:,el2{2}} > 0;
    data = data(ind,:);
    
    yobs = data{:,el2{1}}./data{:,el2{2}};
elseif strcmpi(method,'nocensor')
    ind = data{:,el2} > 0;
    data = data(ind,:);
    
    if strcmp(yscale,'log')
        yobs = log10(data{:,el2});
    else
        yobs = data{:,el2};
    end
else
    yobs = data{:,el2};
end

y = nan(length(x)-nx+1,5);
n = zeros([length(x)-nx+1 1]);

for i = nx/2:length(x)-nx/2
    ind = x(i-nx/2+1) <= data{:,el1} & data{:,el1} <= x(i+nx/2)';
    n(i-nx/2+1) = sum(ind);
    
    switch lower(method)
        case 'nocensor'
            if sum(ind) < 10
               continue;
            end

            [mu,sigma,mu_95,sigma_95] = normfit(yobs(ind),0.05);
            
            [~,~,mu_50,sigma_50] = normfit(yobs(ind),0.50);
                y(i-nx/2+1,:) = [mu_95(1), mu_50(1), mu, mu_50(2), mu_95(2)];
            %switch lower(errtype)
            %    case 'q'
            %        y(i-nx/2+1,:) = mu + sigma*sqrt(2)*erfinv(2*[0.025 0.25 0.5 0.75 0.975] - 1);
            %    case 'qerr'
            %        [~,~,mu_50,sigma_50] = normfit(yobs(ind),0.50);
            %        y(i-nx/2+1,:) = [mu_95(1), mu_50(1), mu, mu_50(2), mu_95(2)];
            %    otherwise
            %        error('Unknown error type.');
            %end
            
        case 'ideal'
            ytmp = prepdata_avg(data(ind,:),el2,yscale);
            y(i-nx/2+1,:) = ytmp(:,1);
        case 'cdf'
            ytmp = prepdata_avg(data(ind,:),el2,yscale);
            y(i-nx/2+1,:) = ytmp(:,2);
        otherwise
            error('Unknown method.');
    end
end

% fix size of x
x = midpt(x);
x = x(nx/2:end-nx/2+1);

if any(strcmpi({'cdf','ideal'},method))
    switch lower(errtype)
        case 'q'
            y = repmat(y(:,3),1,5) + (y - repmat(y(:,3),1,5));
        case 'qerr'
            y = repmat(y(:,3),1,5) + (y - repmat(y(:,3),1,5))./sqrt(repmat(n,1,5));
        otherwise
            error('Unknown error type.');
    end
end

% only needed for histograms
xedges = [x-0.5*(x(2)-x(1)); x(end)+0.5*(x(2)-x(1))];

% remove y values with nan's
ind = all(~isnan(y),2);
x = x(ind);
y = y(ind,:);

% prepare data for plotting
y95 = [y(:,1); flipud(y(:,5)); y(1,5)];
y50 = [y(:,2); flipud(y(:,4)); y(1,2)];
X = [x; flipud(x); x(1)];

% plot histogram?
if plothist == 1
    subplot(3,1,1:2);
    hold on;
end

% uncertainty bounds
if isempty(colour)
    f1 = fill(X,y95, [1 1 1]*.8,'EdgeColor','none');
    hold on
    f2 = fill(X,y50, [1 1 1]*.6,'EdgeColor','none');
else
    f1 = fill(X,y95, colour + (1-colour)*.8,'EdgeColor','none');
    hold on
    f2 = fill(X,y50, colour + (1-colour)*.6,'EdgeColor','none');
end

% transparency
if ~isempty(alpha)
    f1.FaceAlpha = alpha;
    f2.FaceAlpha = alpha;
end

% median line
p = plot(x,y(:,3), 'LineWidth',2);
if ~isempty(colour)
    p.Color = colour;
end

set(gca,'Box','on')
legend({'95% quantile error','50% quantile error','Median'},'Location','northwest')

% ylabel
if strcmp(yscale,'log')
    yl = get(gca,'YLim');
    hpax([floor(yl(1)) ceil(yl(2))],'y');
end
ylabel(ylbl);

% now plot histogram?
if plothist == 1
    subplot(3,1,3);
    hold on;
    histogram(data{:,el1},'BinEdges',xedges,'FaceColor',p.Color,'EdgeColor','none');
    ylabel('No. Analyses');
    set(gca,'Box','on');
end

% x labels
if strcmp(xscale,'log')
    xl = get(gca,'XLim');
    hpax([floor(xl(1)) ceil(xl(2))],'x');
end
xlabel(xlbl);

return