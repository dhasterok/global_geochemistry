function h = scatterlegend(ptsize,ptlabel,varargin)
% scatterlegend - add a point size legend to a scatter plot
%
%   h = scatterlegend(ptsize,ptlabel) will create a legend on a scatterplot
%   with a vector of point sizes, ptsize, with associated labels in a cell
%   array, ptlabel.
%
%   Use the following option value pairs to change the behavior:
%
%       'PointFill'         will fill the points if true, default is false
%
%       'Symbol'            change the symbol using the standard
%                           plot/scatter symbols, default is 'o'
%
%       'Orientation'       change from 'vertical' to 'horizontal'
%                           orientation
%
%       'Location'          default is 'southwest', can also be
%                           'southeast', 'northeast' and 'northwest'
%
%       'Spacing'           change the relative spacing between the points
%                           by supplying a positive value, default is 1.
%                           A spacing > 1 will move the points farther
%                           apart and < 1 closer together

% D. Hasterok, University of Adelaide, 2022

% parse inputs
p = inputParser;
addRequired(p,'ptsize',@isnumeric);
addRequired(p,'ptlabel',@iscell);
addParameter(p,'PointFill',false,@islogical);
addParameter(p,'Symbol','o',@ischar);
addParameter(p,'Orientation','vertical',@ischar);
addParameter(p,'Location','southeast',@ischar);
addParameter(p,'Spacing',1,@isnumeric);

parse(p,ptsize,ptlabel,varargin{:});

np = length(ptsize);
if length(ptlabel) ~= np
    error('The length of points and labels must be the same.');
end

% set properties
sym = p.Results.Symbol;
orientation = p.Results.Orientation;
location = p.Results.Location;
delta = p.Results.Spacing;
ptfill = p.Results.PointFill;

% get current axes limits
xl = get(gca,'XLim');
yl = get(gca,'YLim');

% change orientation and compute initial locations of points and text
switch orientation
    case 'vertical'
        [xp,yp,xt,yt] = ptlocations(np,xl,yl,delta);
        HA = 'left';
        VA = 'middle';
    case 'horizontal'
        % x becomes y and y becomes x
        [yp,xp,yt,xt] = ptlocations(np,xl,yl,delta);
        HA = 'center';
        VA = 'bottom';
    otherwise
        error('Unknown orientation.');
end

% change location and shift points and text if necessary
if contains(location,'east')
    xp = flipud(xl(2) - (xp - xl(1)));
    xt = flipud(xl(2) - (xt - xl(1)));
    if strcmp(orientation,'vertical')
        HA = 'right';
    end
end

if contains(location,'north')
    yp = flipud(yl(2) - (yp - yl(1)));
    yt = flipud(yl(2) - (yt - yl(1)));
    if strcmp(orientation,'horizontal');
        VA = 'top';
    end
end

% plot points
hold on;
if ptfill
    h = scatter(xp,yp,ptsize,sym,'filled');
else
    h = scatter(xp,yp,ptsize,sym);
end

% plot labels
for i = 1:length(ptlabel)
    text(xt(i),yt(i),ptlabel{i},'HorizontalAlignment',HA,'VerticalAlignment',VA);
end

return


function [xp,yp,xt,yt] = ptlocations(np,xl,yl,delta)

dx = diff(xl);
dy = diff(yl);

xp = repmat(xl(1) + 0.05*dx,np,1);
yp = linspace(yl(1) + 0.05*dy,delta*(yl(1)+0.3*dy),np)';

xt = xp + 0.02*dx;
yt = yp;

return