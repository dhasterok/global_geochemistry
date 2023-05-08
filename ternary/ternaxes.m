function ternaxes(varargin)
% TERNAXES - Creates ternary axes
%
%   ternaxes will create a basic ternary axes
%
%   ternaxes('NVertices',n) where n = 3 will produce a single ternary axis,
%   and n = 4 will create a double (mirrored) set of ternary axes.
%
%   n = 3           n = 4
%           A               A
%          / \             / \
%        B --- C         B --- C
%                          \ /
%                           D
%
%   Options:
%
%       'NVertices'             Number of vertices, 3 or 4
%
%       'Axes'                  Axes object for plotting, default = gca
%
% See also ternary, terngrid, terntick

% Last Modified: 8 May 2023
% D. Hasterok, University of Adelaide

% parse inputs
% ------------------------
p = inputParser;

addParameter(p,'NVertices',3,@isnumeric);   % number of axes
addParameter(p,'Axes',[],@isgraphics);      % handle to axes

parse(p,varargin{:});

n = p.Results.NVertices;
if n ~= 3 && n ~= 4
    error('NVertices can be 3 or 4.');
end

ax = p.Results.Axes;
if isempty(ax)
    ax = gca;
end
% ------------------------

w = 0.5;            % half width
h = 0.5/tan(pi/6);  % vertical scale

hold(ax,'on');

%create axes
fill(ax,[-w 0 w -w],[0 h 0 0],'w');
plot(ax,[-w 0 w -w],[0 h 0 0],'k-','LineWidth',1);
if n == 4
    fill(ax,[-w 0 w -w],-[0 h 0 0],'w');
    plot(ax,[-w 0 w -w],-[0 h 0 0],'k-','LineWidth',1);
end

ax.Visible = 'off';
ax.DataAspectRatio = [1 1 1];
axis(ax,'tight');

return