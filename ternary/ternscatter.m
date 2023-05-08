function [t,varargout] = ternscatter(varargin)
% TERNSCATTER - plot data on a ternary axes
%
%    t = ternplot(a,b,c) or
%    t = ternplot(a,b,c, sym) where sym is the line specification, see
%    PLOT.
%
%    To use ternplot like scatter (points colored using a fourth vector):
%    t = ternplot(a,b,c, sym, {S,P}) where S is the size and P is the 
%    vector of values to be scaled as the color.  S and P can be a single
%    value or of length equal to the number of points.  S and P lengths
%    need not be the same.
%
%    The order of the axes are as follows:
%
%                    A
%                   / \
%                  /   \
%                 B --- C
%
%    See also ternary, ternplot

% Last Modified: 8 May 2023
% D. Hasterok, University of Adelaide

% parse inputs
% ------------------------
p = inputParser;

% data vectors for plotting on ternary
addRequired(p,'A',@isnumeric);
addRequired(p,'B',@isnumeric);
addRequired(p,'C',@isnumeric);
addOptional(p,'D',@isnumeric);

% point properties
addParameter(p,'Size',36,@isnumeric);           % size
addParameter(p,'Color',lines(1),@isnumeric);    % color
addParameter(p,'Group',{});                     % group (categorized)
addParameter(p,'Alpha',1,@isnumeric);           % transparency
addParameter(p,'Symbol','o',@ischar);           % symbol

addParameter(p,'Axes',[],@isgraphics);          % handle to axes
parse(p,varargin{:});

a = p.Results.A;
b = p.Results.B;
c = p.Results.C;
d = p.Results.D;

ax = p.Results.Axes;
if isempty(ax)
    ax = gca;
end

% point properties
S = p.Results.Size;
C = p.Results.Color;
G = p.Results.Group;

% transparency
alpha = p.Results.Alpha;

% symbols
sym = p.Results.Symbol;
% ------------------------


if isempty(d)
    [x,y] = tern2xy(a,b,c);
else
    x = zeros(size(a));
    y = zeros(size(a));
    
    ind = d > 0;
    [x(~ind),y(~ind)] = tern2xy(a(~ind),b(~ind),c(~ind));
    [x(ind),y(ind)] = tern2xy(d(ind),b(ind),c(ind));
    y(ind) = -y(ind);
end  

if isempty(G)
    t = scatter(ax, x,y, S,C, sym,'filled', 'MarkerFaceAlpha',alpha);
else
    %group = unique(C);
    %t = gscatter(ax, x,y, G,C, sym, sqrt(S));
    [t,leg] = catscatter(x,y,G, 'Size',S, 'Symbol',sym, 'Alpha',alpha, 'Axes',ax);
    varargout{2} = leg;
end

return
