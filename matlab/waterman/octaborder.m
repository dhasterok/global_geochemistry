function varargout = octaborder(varargin)
% OCTABORDER - produces a border for an octahedral plot
%
%   h = octaborder(dl) produces a border for an octahedral map.  Graphics
%   handles for the graticules are returned in the array h.
%
%   There are several option value pairs:
%
%       'Projection'        choose between Waterman and Cahill
%                           Cahill not yet implemented
%
%       'Style'             'Butterfly' or 'M'
%
%       'Shift'             shifts the map latitude to the east by the 
%                           specified angle in degrees
%
%       'Inset'             plot Antarctica as contiguous, true or false
%                           default is false
%
%       'MaxLat'            if an inset is included, then a maximum
%                           latitude extent for the inset may be included.
%                           default is -60, but -50 will capture the entire
%                           Scotia Plate
%
%       'Export'            if true, will export as a shapefile
%
% see also OCTAHEDRALPROJECTION, OCTAGRATICULE, OCTACOAST, OCTARASTER

% parse inputs
p = inputParser;
addParameter(p,'Projection','Waterman',@ischar);
addParameter(p,'Style','Butterfly',@ischar);
addParameter(p,'Inset',false,@islogical);
addParameter(p,'MaxLat',-60,@isnumeric);
addParameter(p,'Export',false,@islogical);
addParameter(p,'Scale',1,@isnumeric);

parse(p,varargin{:});
scale = p.Results.Scale;

ax = gca;
hold(ax,'on');
prepmap;

%xmin = inf;
%xmax = -inf;
%ymin = inf;
%ymax = -inf;

% plot border
if p.Results.Inset
    % with Antarctic inset
    [xb,yb,xa,ya] = panelborders(p);
    xb = scale*xb;
    yb = scale*yb;
    xa = scale*xa;
    ya = scale*ya;

    h = plot(ax,xb,yb,'k-','LineWidth',1);
    h(2) = plot(ax,xa,ya,'k-','LineWidth',1);
else
    % without inset
    [xb,yb] = panelborders(p);
    xb = scale*xb;
    yb = scale*yb;
    
    h = plot(ax,xb,yb,'k-','LineWidth',1);
end

if nargout == 1
    varargout{1} = h;
end

% if Export is false, return
if ~p.Results.Export
    return
end

% export border to shapefile
filename = ['border_',p.Results.Projection,'_',p.Results.Style];

xb(end+1,:) = NaN;
yb(end+1,:) = NaN;

if p.Results.Inset
    filename = [filename,'_Inset.shp'];
    x = [xb(:); xa];
    y = [yb(:); ya];
else
    filename = [filename,'.shp'];
    x = xb(:);
    y = yb(:);
end

s = mapshape(x,y);
s.Geometry = 'line';

shapewrite(s,filename);

return