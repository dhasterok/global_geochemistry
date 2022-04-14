function varargout = cloudcity(val,varargin)
% CLOUDCITY - Creates mirrored histograms (inspired by Star Wars)
%
%   [h,edges,n] = cloudcity(x,N) plots a histogram with bars oriented
%   horizontally rather than vertically and mirrored.  Return values
%   include h, the graphics handle, the bin edges, and n, the data in each
%   bin.
%
%
%   Input additional options as cloudcity(data, 'option', value, ...).
%       BinEdges    Edges of the histogram bins, default creates
%                   sqrt(length(x)) bins starting at the lowest value of x
%                   and ending at the largest value of x
%   	Scale       Normalization factor for plotting, default normalizes
%                   such that largest bin is 1
%       Shift       Shift applied to the zero value, useful for plotting
%                   mulitiple histograms
%       RefAxis     Axis of reflection for the mirroring, default is 'y'

% Last updated: 24 July 2020

if isempty(val)
    warning('Value vector was empty.');
    return;
end

% parse inputs
shift = 0;
coord = 'y';
edges = linspace(min(val),max(val),round(sqrt(length(val))));
N = 1;
if nargin > 2
    opt = 1;
    while opt + 1 < nargin
        switch lower(varargin{opt})
            case 'scale'
                N = varargin{opt+1};
            case 'shift'
                shift = varargin{opt+1};
            case 'refaxis'
                coord = varargin{opt+1};
            case 'binedges'
                edges = varargin{opt+1};
            otherwise
                error('Unknown option.');
        end
        opt = opt + 2;
    end
end

% create x values for plotting bins
x = [edges; edges];
x = x(:);
x = [x; flipud(x)];

% compute number in each bin
n = histc(val,edges)';

% create y values for plotting
y = [0 n(1:end-1); n(1:end-1) 0]*N/max(n);
y = y(:);

% shift data
y = [y; -flipud(y)] + shift;

% axis for plotting mirror (swap x,y) if necessary
if strcmp(coord,'y');
    h = fill(y,x,[0.7 0.7 0.7]);
else
    h = fill(x,y,[0.7 0.7 0.7]);
end

% prepare outputs
if nargout >= 1
    varargout{1} = h;
    if nargout <= 3 
        varargout{2} = edges;
        varargout{3} = n;
    elseif nargout > 3
        error('Too many outputs requested');
    end
end

return
