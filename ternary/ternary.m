function ternary(varargin)
% TERNARY - creates ternary axes.
%
%   ternary will create a basic set of ternary axes with a grid and tick
%   marks.
%
%   To add labels, ternary('Labels',{'A','B','C'}) will additionally label
%   the ternary axes.  Add an additional label to create a diamond,
%   ternary('Labels',{'A','B','C','D'}).
%
%      triangular      diamond
%
%           A              A
%          / \            / \
%        B --- C        B --- C
%                         \ /
%                          D
%
%   Options:
%
%       'Labels'                Axes labels as a cell array of
%                               strings/chars
%
%       'GridSpacing'           Spacing of grid lines, 0 = no grid,
%                               default = 0.1
%
%       'TickSpacing'           Spacing of tick lines, 0 = no tick,
%                               default = 0.1
%
%       'Axes'                  Axes object for plotting, default = gca
%
% See also ternaxes, terngrid, terntick

% Last Modified: 8 May 2023
% D. Hasterok, University of Adelaide

% parse inputs
% ------------------------
p = inputParser;

addParameter(p,'Labels',{},@iscellstr);            % labels
addParameter(p,'GridSpacing',0.1,@isnumeric);   % grid spacing
addParameter(p,'TickSpacing',0.1,@isnumeric);   % grid spacing
addParameter(p,'Axes',[],@isgraphics);          % handle to axes

parse(p,varargin{:});

axlbl = p.Results.Labels;
dg = p.Results.GridSpacing;
dt = p.Results.TickSpacing;

ax = p.Results.Axes;
if isempty(ax)
    ax = gca;
end
hold(ax,'on');
% ------------------------

% create ternary axes
len = length(axlbl);
if len == 0
    len = 3;
end
ternaxes('NVertices',len,'Axes',ax);

% add grid lines
if dg ~= 0
    terngrid('NVertices',len,'Axes',ax,'Spacing',dg);
end

% add tick marks
if dt ~= 0
    terntick('NVertices',len,'Axes',ax,'Spacing',dt);
end

% add labels
if length(axlbl) == 3
    axis(ax, [-1.346255486272525 1.346255486272525 -0.480230082488086 1.346255486272525]);
    ternlabel(axlbl{1},axlbl{2},axlbl{3}, 'Axes',ax);
elseif length(axlbl) == 4
    axis(ax, [-1.346255486272525 1.346255486272525 -1.346255486272525 1.346255486272525]);
    ternlabel(axlbl{1},axlbl{2},axlbl{3},axlbl{4}, 'Axes',ax);
end

ax.DataAspectRatio = [1 1 1];
axis(ax,'tight');

return