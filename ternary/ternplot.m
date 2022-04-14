function t = ternplot(a,b,c,varargin)
% TERNPLOT - plot data on a ternary axes
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
%    See also TERNARY

opt = 1;
var = 3;
if isnumeric(varargin{1});
    d = varargin{1};
    opt = 2;
    var = 4;
end
switch nargin - opt + 1
    case 5
        sym = varargin{opt};
        lspec = varargin{opt+1};
        flag = 1;
    case 4
        sym = varargin{opt};
        [nr,nc] = size(sym);
        if ischar(sym) & nr == 1
            flag = 0;
        elseif nr > 1
            flag = 1;
            lspec = {ones(size(a)), 'o'};
        else
            error(['Fourth argument must be a string with a single ', ...
                'row LINESPEC, or a cell of LINESPEC with the same ', ...
                'dimensions as input points.']);
        end
    case 3
        sym = [];
    otherwise
        error('Incorrect number of inputs.');
end

if var == 3
    [x,y] = tern2xy(a,b,c);
elseif var == 4
    x = zeros(size(a));
    y = zeros(size(a));
    
    ind = d > 0;
    [x(~ind),y(~ind)] = tern2xy(a(~ind),b(~ind),c(~ind));
    [x(ind),y(ind)] = tern2xy(d(ind),b(ind),c(ind));
    y(ind) = -y(ind);
end  

if ~isempty(sym)
    if flag
        t = scatter(x,y,lspec{1},lspec{2},sym,'filled');
        %set(t,'MarkerEdgeColor','k');
    else
        t = plot(x,y,sym);
    end
else
    t = plot(x,y);
end

return
