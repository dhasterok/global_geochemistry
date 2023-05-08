function ternlabel(varargin)
% ternlabel - adds labels to ternary axes.
%
%   ternlabel(A,B,C) will label the ternary axes.
%
%           A
%          / \
%        B --- C
%
%   ternlabel(A,B,C,D) will label a double ternary axes.
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

% D. Hasterok, University of Adelaide
axlbl = {};
ax = [];
opt = 1;
c = 1;
while opt < nargin + 1
    switch varargin{opt}
        case 'Axes'
            ax = varargin{opt+1};
            if ~isgraphics(ax,'Axes');
                error('Argument must be an axes object.');
            end
            opt = opt + 2;
        otherwise
            if c > 4
                error('Too many axes inputs.');
            elseif ischar(varargin{opt}) | isstring(varargin{opt})
                axlbl{c} = varargin{opt};
            else
                error('Expected label should be of type char or string.')
            end
            c = c + 1;
            opt = opt + 1;
    end
end
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
