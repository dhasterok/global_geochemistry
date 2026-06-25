function varargout = octapolyplot(tmp,varargin)

fileflag = false;

% parse inputs
p = inputParser;

if isnumeric(tmp)
    lat = tmp;
    clear tmp

    addRequired(p,'Longitude',@numeric);
else
    filename = tmp;
    clear tmp;

    switch lower(filename(end-3:end))
        case '.shp'
            s = shaperead(filename);
        otherwise
            error('Unknown filetype.');
    end
    fileflag = true;
end
    
addParameter(p,'Projection','Waterman',@ischar);
addParameter(p,'Style','Butterfly',@ischar);
addParameter(p,'Shift',0,@isnumeric);
addParameter(p,'Transform','forward',@ischar);
addParameter(p,'Inset',false,@islogical);
addParameter(p,'MaxLat',-50,@isnumeric);
addParameter(p,'Export',false,@islogical);

parse(p,varargin{:});

if fileflag
    if p.Results.Inset
        sa = s;
        ind = true(size(s));
        parfor i = 1:length(s)
            [s(i).X,s(i).Y,sa(i).X,sa(i).Y] = trimpoly(s(i).Y,s(i).X,p);
            if isempty(sa(i).X)
                ind(i) = false;
            end
        end
        sa = sa(ind);
    else
        parfor i = 1:length(s)
            [s(i).X,s(i).Y,] = trimpoly(s(i).Y,s(i).X,p);
        end
    end
else
    lon = p.Results.Longitude;
    [x,y] = trimpoly(lat,lon,p.Results.Inset);
end

if p.Results.Export
    filename = [filename(1:end-4),'_Projection_',p.Results.Projection,'_Style_',p.Results.Style,'_Shift_',num2str(p.Results.Shift)];
    shapewrite(s,[filename,'.shp']);
    if p.Results.Inset
        shapewrite(sa,[filename,'_Inset.shp']);
    end
end

return

function [X,Y,varargout] = trimpoly(lat,lon,p)

persistent quadpoly insetpoly

if ~exist('quadpoly','var') | isempty(quadpoly)
    % create polygons for clipping
    if p.Results.Inset
        [xb,yb,xbi,ybi] = panelborders(p);
    
        insetpoly = polyshape(xbi,ybi,'Simplify',true);
    else
        [xb,yb] = panelborders(p);
    end

    quadpoly = [];
    for n = 1:4
        quadpoly = [quadpoly; polyshape(xb(:,n),yb(:,n),'Simplify',true)];
    end
end

ax = gca;
hold(ax,'on');
prepmap;

if isempty(lat)
    X = [];
    Y = [];
    return;
end
    
% split lines at map octahedral boundaries
[x,y,xa,ya] = OctahedralProjection(lat,lon, ...
    'Projection',p.Results.Projection, ...
    'Style',p.Results.Style, ...
    'Shift',p.Results.Shift, ...
    'Inset',p.Results.Inset, ...
    'MaxLat',p.Results.MaxLat);

x = [NaN; x'];
y = [NaN; y'];

X = [];
Y = [];
ind = find(isnan(x));
for k = 2:length(ind)
    for n = 1:4
        poly = polyshape([x(ind(k-1)+1:ind(k)-1); x(ind(k-1)+1)], [y(ind(k-1)+1:ind(k)-1); y(ind(k-1)+1)],'Simplify',true);

        % find intersection of polygons
        polyi = intersect(poly,quadpoly(n));

        if ~isempty(polyi.Vertices)
            plot(ax,polyi,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.4,'LineWidth',0.25,'EdgeColor','none');
            
            X = [X; NaN; polyi.Vertices(:,1)];
            Y = [Y; NaN; polyi.Vertices(:,2)];
        end
    end
end

if p.Results.Inset% && ~isempty(xa)
    xa = [NaN; xa];
    ya = [NaN; ya];
    
    Xa = [];
    Ya = [];
    ind = find(isnan(xa));
    for k = 2:length(ind)
        poly = polyshape([xa(ind(k-1)+1:ind(k)-1); xa(ind(k-1)+1)], [ya(ind(k-1)+1:ind(k)-1); ya(ind(k-1)+1)]);

        % find intersection of polygons
        polyi = intersect(poly,insetpoly);

        if ~isempty(polyi.Vertices)
            plot(ax,polyi,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.4,'LineWidth',0.25,'EdgeColor','none');
            
            Xa = [Xa; NaN; polyi.Vertices(:,1)];
            Ya = [Ya; NaN; polyi.Vertices(:,2)];
        end
    end
    varargout{1} = Xa;
    varargout{2} = Ya;
end

return


