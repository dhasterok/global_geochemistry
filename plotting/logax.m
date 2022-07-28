function logax(lim,varargin)
% logax - Produces log-axes limits
%
%   logax(lim) where lim are the log10 values of the axes limits will
%   produce log-axes tick marks with powers of 10 values in decimal form.
%   Yes, I know matlab can do this in other ways, but this function
%   produces labels the tick marks the way I like them.  It will also add
%   them to plots that may not ordinarly allow changing from linear to log
%   axes.
%
%   The following option-value pairs may be given:
%
%       Axis                'x' or 'y' to add ticks to x- or y-axis,
%                           default is y
%
%       Label               add an axis label
%
%       TickLabelRotation   angle of text rotation, default is 0

% Last modified 6 July, 2022
% D. Hasterok, University of Adelaide

p = inputParser;
addRequired(p,'lim',@isnumeric);
addParameter(p,'Axis','y',@ischar);
addParameter(p,'Label','',@ischar);
addParameter(p,'TickLabelRotation',0,@isnumeric);

parse(p,lim,varargin{:});
lim = p.Results.lim;
x_or_y = p.Results.Axis;
axlbl = p.Results.Label;
angle = p.Results.TickLabelRotation;

% create cell array of text
mt = log10([1:9]);
tmp = {'','','','','','','',''};
atic = [];
aticlbl = {};
for i = lim(1):lim(2)
    atic = [atic,i + mt];
    aticlbl = [aticlbl,{num2str(10^i)},tmp];
end

switch lower(x_or_y)
    case 'x'
        % add optional label
        if ~isempty(axlbl)
            xlabel(axlbl);
        end
        %set(gca,'XTick',[lim(1):lim(2)],'XTickLabel',10.^[lim(1):lim(2)],'Box','on');
        set(gca,'XTick',atic,'XTickLabel',aticlbl,'Box','on','XTickLabelRotation',angle);
        xlim(lim);
    case 'y'
        % add optional label
        if ~isempty(axlbl)
            ylabel(axlbl);
        end
        %set(gca,'YTick',[lim(1):lim(2)],'YTickLabel',10.^[lim(1):lim(2)],'Box','on');
        set(gca,'YTick',atic,'YTickLabel',aticlbl,'Box','on','YTickLabelRotation',angle);
        ylim(lim);
    otherwise
        warning('Ignoring hpax command, incorrect axes argument');
end

return