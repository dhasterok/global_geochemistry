function varargout = octacoast(varargin)
% octacoast - produces a coastline on an octahedral plot
%
%   h = octacoast(dl) produces a coastline on an octahedral map.  A
%   graphics handles for the coastline is returned in h.
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
% see also OCTAHEDRALPROJECTION, OCTAGRATCULE, OCTABORDER, OCTARASTER

addpath ../maptools

% parse inputs, they are passed straight through to OctahedralProjection
p = inputParser;
addParameter(p,'Projection','Waterman',@ischar);
addParameter(p,'Style','Butterfly',@ischar);
addParameter(p,'Shift',0,@isnumeric);
addParameter(p,'Inset',false,@islogical);
addParameter(p,'MaxLat',-60,@isnumeric);
addParameter(p,'Export',false,@islogical);
addParameter(p,'Scale',1,@isnumeric)

parse(p,varargin{:});

scale = p.Results.Scale;
% found in maptools/
load('coast.mat');

[x,y,xa,ya] = cutpolyconvert(clat,clon,p);
x = scale*x;
y = scale*y;
xa = scale*xa;
ya = scale*ya;

% plot coastline
ax = gca;
hold(ax,'on');
prepmap;
h = plot(ax,x,y,'k-');

if p.Results.Inset
    plot(ax,xa,ya,'k-');
end

if nargout == 1
    varargout{1} = h;
end

if ~p.Results.Export
    return
end

filename = ['coast_',p.Results.Projection,'_',p.Results.Style,'_Shift_',num2str(p.Results.Shift)];
if p.Results.Inset
    filename = [filename,'_Inset.shp'];
    x = [x; NaN; xa];
    y = [y; NaN; ya];
else
    filename = [filename,'.shp'];
end

s = mapshape(x,y);
s.Geometry = 'line';

shapewrite(s,filename);

return