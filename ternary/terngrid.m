function terngrid(varargin)
% TERNGRID - Adds ternary grid lines to ternary axes
%
%   terngrid will create a basic ternary grid every 0.1 units.
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
% See also ternary, ternaxes, terntick

% Last Modified: 8 May 2023
% D. Hasterok, University of Adelaide

% parse inputs
% ------------------------
p = inputParser;

addParameter(p,'Spacing',0.1,@isnumeric);   % grid spacing
addParameter(p,'NVertices',3,@isnumeric);   % number of axes
addParameter(p,'Axes',[],@isgraphics);      % handle to axes

parse(p,varargin{:});

dg = p.Results.Spacing;

n = p.Results.NVertices;
if n ~= 3 && n ~= 4
    error('NVertices can be 3 or 4.');
end

ax = p.Results.Axes;
if isempty(ax)
    ax = gca;
end
% ------------------------

% grid
xa = [-0.5+dg:dg:0.5-dg];
ya = zeros(size(xa));

[a,b,c] = xy2tern(xa,ya);
[xb,yb] = tern2xy(b,a,c);
[xc,yc] = tern2xy(c,b,a);

plot(ax, [xa;xb],[ya;yb],'-','LineWidth',0.25,'Color',[0.8 0.8 0.8]);
plot(ax, [xb;fliplr(xc)],[yb;fliplr(yc)],'-','LineWidth',0.25,'Color',[0.8 0.8 0.8]);
plot(ax, [xc;xa],[yc;ya],'-','LineWidth',0.25,'Color',[0.8 0.8 0.8]);

if n == 4
    plot(ax, [xa;xb],-[ya;yb],'-','LineWidth',0.25,'Color',[0.8 0.8 0.8]);
    plot(ax, [xb;fliplr(xc)],-[yb;fliplr(yc)],'-','LineWidth',0.25,'Color',[0.8 0.8 0.8]);
    plot(ax, [xc;xa],-[yc;ya],'-','LineWidth',0.25,'Color',[0.8 0.8 0.8]);
end

return

