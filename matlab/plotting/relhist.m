function h = relhist(x,y,varargin)


Cleft = [227, 179, 179]/255;
Cright = [159, 209, 153]/255;

xl = [];
bw = [];
edges = [];
hinge = 0;
normtype = '';

opt = 1;
while opt + 2 <= nargin
    switch lower(varargin{opt})
        case 'binwidth'
            bw = varargin{opt+1};
            opt = opt + 2;
        case 'binedges'
            edges = varargin{opt+1};
            opt = opt + 2;
        case 'normalization'
            normtype = varargin{opt+1};
            opt = opt + 1;
        case 'leftcolor'
            Cleft = varargin{opt+1};
            opt = opt + 2;
        case 'rightcolor'
            Cright = varargin{opt+1};
            opt = opt + 2;
        case 'hinge'
            hinge = varargin{opt+1};
            opt = opt + 2;
        case 'xlim'
            xl = varargin{opt+1};
            opt = opt + 2;
        otherwise
            error('Unknown option.');
    end
end

r = x - y;

rind = r > 0 + hinge;
lind = r <= 0 + hinge;

ax = gca;
hold(ax,'on');
ax.Box = 'on';

h(1) = histogram(ax,r(lind));
h(1).FaceColor = Cleft;
h(1).EdgeColor = 'none';
h(2) = histogram(ax,r(rind));
h(2).FaceColor = Cright;
h(2).EdgeColor = 'none';

if ~isempty(bw)
    h(1).BinWidth = bw;
    h(2).BinWidth = bw;
elseif ~isempty(edges)
    h(1).BinEdges = edges;
    h(2).BinEdges = edges;
end

if ~isempty(normtype)
    h(1).Normalization = normtype;
    h(2).Normalization = normtype;
end

if ~isempty(xl)
    ax.XLim = xl;
end

yl = ax.YLim;
h(3) = plot(ax,[0 0],yl,'k-','LineWidth',1);

return