function ternlabel(varargin)
% ternlabel - adds labels to ternary axes.
%
%   ternlabel({'A','B','C'}) will label the ternary axes from a cell array.
%
%           A
%          / \
%        B --- C
%
%   ternlabel({'A','B','C','D'}) will label a double ternary axes.
%
%           A
%          / \
%        B --- C
%          \ /
%           D
%
%   Options:
%
%       'Axes'                  Axes object for plotting, default = gca
%
% See also ternaxes, terngrid

% Last Modified: 8 May 2023
% D. Hasterok, University of Adelaide

% parse inputs
% ------------------------
p = inputParser;

addRequired(p,'Labels',@iscellstr);     % labels
addParameter(p,'Axes',[],@isgraphics);      % handle to axes

parse(p,varargin{:});

axlbl = p.Results.Labels;
n = length(axlbl);

ax = p.Results.Axes;
if isempty(ax)
    ax = gca;
end
% ------------------------

w = 0.5;
h = 0.5/tan(pi/6);
d = 0.02;

text(ax, 0,h+d,axlbl{1},'FontSize',12,'HorizontalAlignment','center','FontWeight','bold','VerticalAlignment','bottom');
if length(axlbl) == 3
    text(ax, -w,-d,axlbl{2},'FontSize',12,'HorizontalAlignment','center','FontWeight','bold');
    text(ax, w,-d,axlbl{3},'FontSize',12,'HorizontalAlignment','center','FontWeight','bold');
else
    text(ax, -w-d,0,axlbl{2},'FontSize',12,'HorizontalAlignment','right','FontWeight','bold');
    text(ax, w+d,0,axlbl{3},'FontSize',12,'HorizontalAlignment','left','FontWeight','bold');
    text(ax, 0,-(h+d),axlbl{4},'FontSize',12,'HorizontalAlignment','center','FontWeight','bold','VerticalAlignment','top');
end

return
