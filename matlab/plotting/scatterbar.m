function varargout = scatterbar(x, y, varargin)
% SCATTERBAR - Makes error bar plots using scatter.
%   h = scatterbar(x,y,...,ebar) plots data over error bars.  Error bars
%   (ebar) can be specified for 'x', 'y', or 'both'.  An ebar option must
%   be provided.
%
%   As with scatter plots, points can be sized and colored by including
%   'Size' and 'Color' options.  Options are supplied as option-value
%   pairs, e.g., scatterbar(x,y,...,ebar,'Option',value,...).
%
%   Error bars can be symmetric or asymmetric about the value by supplying
%   a single bar value or two bar values for each point.  For example to
%   create asymmetric y-error bars, scatterbar(x,y,ey_lower,ey_upper,'y').
%
%   Input options:
%       Size        Size of markers.
%       Color       Color of markers.
%       Shape       Marker shape (same as plot).
%       Type        Type of data supplied for error bars, either 'sigma' as
%                   in x +/- sigma or y +/- sigma.  To supply or the 
%                   coordinates of the end points, 'value' for
%                   (xleft,xright) or (ylower,yupper).

% convert all data to rows (most convenient for plotting errorbars)
x = row(x);
y = row(y);

% find the error bar type (x, y or both)
for i = 2:nargin-2
    if any(strcmp(class(varargin{i}),{'char','string'}))
        ebar = varargin{i};
        c = i+2;
        break;
    end
end

% read read error bar data
switch ebar
    case 'x'
        ey = [];
        if c == 4
            ex = [row(varargin{1}); row(varargin{1})];
        elseif c == 5
            ex = [row(varargin{1}); row(varargin{2})];
        else
            error('Too many errorbars.');
        end
    case 'y'
        ex = [];
        if c == 4
            ey = [row(varargin{1}); row(varargin{1})];
        elseif c == 5
            ey = [row(varargin{1}); row(varargin{2})];
        else
            error('Too many errorbars.');
        end
    case 'both'
        if c == 5
            ex = [row(varargin{1}); row(varargin{1})];
            ey = [row(varargin{2}); row(varargin{2})];
        elseif c == 7
            ex = [row(varargin{1}); row(varargin{2})];
            ey = [row(varargin{3}); row(varargin{4})];
        else
            error('Incorrect number of errorbars.');
        end
    otherwise
        error('Incorrect error bar coordinates.');
end

% default options
marker_size = [];
marker_color = [];
marker_shape = [];
bar_type = 'sigma';

% read options
opt = c + 1;
while opt < nargin
    switch lower(varargin{opt})
        case 'size'
            marker_size = varargin{opt+1};
        case 'color'
            marker_color = varargin{opt+1};
        case 'shape'
            marker_shape = varargin{opt+1};
        case 'type'
            bar_type = varargin{opt+1};
    end
    opt = opt + 2;
end

% make errorbars for plotting
xbar = ebar_type(x,ex,bar_type);
ybar = ebar_type(y,ey,bar_type);

% store hold state and turn on if necessary
hold_state = ishold;
if ~hold_state
    hold on;
end

% plot x error bars
if ~isempty(xbar)
    h{3} = plot(xbar,[y; y],'-','Color',[0.8,0.8,0.8]);
    %set(get(get(h{3},'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end


% plot y error bars
if ~isempty(ybar)
    h{2} = plot([x; x],ybar,'-','Color',[0.8,0.8,0.8]);
    %set(get(get(h{2},'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end

% make scatter plot (add size and color if necessary)
if isempty(marker_color)
    h{1} = scatter(x',y',marker_size,'filled');
else
    h{1} = scatter(x',y',marker_size,marker_color,'filled');
end

% Change marker from 'o' to marker_shape
if ~isempty(marker_shape)
    h{1}.Marker = marker_shape;
end

% if graphics handles are requested
if nargout == 1
    varargout{1} = h;
end

% if hold was off, turn off again
if ~hold_state
    hold off;
end

return

% prepares error bars for plotting
function ebar = ebar_type(v,ev,type)

if isempty(ev)
    ebar = [];
    return
end

switch lower(type)
    case 'sigma'
        ebar = [v; v] + [-ev(1,:); ev(2,:)];
    case 'value'
        ebar = ev;
    otherwise
        error('Unknown error bar type.');
end

return

% ensure vector is in row form.
function v = row(v)

v = v(:)';

return
