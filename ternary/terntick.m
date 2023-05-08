function terntick(varargin)
% TERNTICK - Adds ternary tick marks to ternary axes
%
%   terntick will create a basic ternary tick marks every 0.1 units.
%
%   Options:
%
%       'Spacing'               Grid spacing, default is 0.1 units
%
%       'NVertices'             Number of vertices, 3 or 4
%
%                               n = 3           n = 4
%                                          A               A
%                                         / \             / \
%                                       B --- C         B --- C
%                                         \ /
%                                          D
%
%       'Axes'                  Axes object for plotting, default = gca
%
% See also ternary, ternaxes, terngrid

% Last Modified: 8 May 2023
% D. Hasterok, University of Adelaide

% parse inputs
% ------------------------
p = inputParser;

addParameter(p,'Spacing',0.1,@isnumeric);   % grid spacing
addParameter(p,'NVertices',3,@isnumeric);   % number of axes
addParameter(p,'Axes',[],@isgraphics);      % handle to axes

parse(p,varargin{:});

dt = p.Results.Spacing;

n = p.Results.NVertices;
if n ~= 3 && n ~= 4
    error('NVertices can be 3 or 4.');
end

ax = p.Results.Axes;
if isempty(ax)
    ax = gca;
end
% ------------------------

% tick marks
xp = [-0.5+dt:dt:0.5-dt];
ya = zeros([3,length(xp)]);
ya([1 3],:) = 0.015;
xa = zeros([3,length(xp)]);
xa(2,:) = xp;
xa(1,:) = xa(2,:) + ya(1,:)/tan(pi/3);
xa(3,:) = xa(2,:) - ya(1,:)/tan(pi/3);

plot(ax, xa,ya,'k-','LineWidth',1);
[a,b,c] = xy2tern(xa,ya);
[xb,yb] = tern2xy(b,a,c);
plot(ax, xb,yb,'k-','LineWidth',1);
[xc,yc] = tern2xy(c,b,a);
plot(ax, xc,yc,'k-','LineWidth',1);

if n == 4
    plot(ax, xa,-ya,'k-','LineWidth',1);
    plot(ax, xb,-yb,'k-','LineWidth',1);
    plot(ax, xc,-yc,'k-','LineWidth',1);
end

return
