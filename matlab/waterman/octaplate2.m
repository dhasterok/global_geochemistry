function varargout = octaplatebndry(varargin)

% parse inputs
p = inputParser;
addParameter(p,'Projection','Waterman',@ischar);
addParameter(p,'Style','Butterfly',@ischar);
addParameter(p,'Shift',0,@isnumeric);
addParameter(p,'Transform','forward',@ischar);
addParameter(p,'Inset',false,@islogical);
addParameter(p,'Dataset','plate',@ischar);

parse(p,varargin{:});

% create polygons for clipping
blon = [-179.999999  -89.999999  0.000001  90.000001;
    -179.999999  -89.999999  0.000001  90.000001;
    -179.999999  -89.999999  0.000001  90.000001;
    -179.999999  -89.999999  0.000001  90.000001;
    -161.3603887  -71.36038869  18.63961131  108.6396113;
    -179.999999  -89.999999  0.000001  90.000001;
    -179.999999  -89.999999  0.000001  90.000001;
    -179.999999  -89.999999  0.000001  90.000001;
    -179.999999  -89.999999  0.000001  90.000001;
    -90.000001  -0.000001  89.999999  179.999999;
    -90.000001  -0.000001  89.999999  179.999999;
    -90.000001  -0.000001  89.999999  179.999999;
    -90.000001  -0.000001  89.999999  179.999999;
    -108.6396113  -18.63961131  71.36038869  161.3603887;
    -90.000001  -0.000001  89.999999  179.999999;
    -90.000001  -0.000001  89.999999  179.999999;
    -90.000001  -0.000001  89.999999  179.999999;
    -90.000001  -0.000001  89.999999  179.999999];

blat = [89.999999  89.999999  89.999999  89.999999;
    71.38865281  71.38865281  71.38865281  71.38865281;
    18.7478136  18.7478136  18.7478136  18.7478136;
    0.000001  0.000001  0.000001  0.000001;
    0  0  0  0;
    -0.000001  -0.000001  -0.000001  -0.000001;
    -18.7478136  -18.7478136  -18.7478136  -18.7478136;
    -71.38865281  -71.38865281  -71.38865281  -71.38865281;
    -89.999999  -89.999999  -89.999999  -89.999999;
    -89.999999  -89.999999  -89.999999  -89.999999;
    -71.38865281  -71.38865281  -71.38865281  -71.38865281;
    -18.7478136  -18.7478136  -18.7478136  -18.7478136;
    -0.000001  -0.000001  -0.000001  -0.000001;
    0  0  0  0;
    0.000001  0.000001  0.000001  0.000001;
    18.7478136  18.7478136  18.7478136  18.7478136;
    71.38865281  71.38865281  71.38865281  71.38865281;
    89.999999  89.999999  89.999999  89.999999];

if p.Results.Inset
    [xb,yb,xa,ya] = panelborders(p);

    polyinset = polyshape(xa,ya,'Simplify',true);
else
    [xb,yb] = panelborders(p);
end

for n = 1:4
    quadpoly(n) = polyshape(blon(:,n),blat(:,n),'Simplify',true);
end

% read polygon shapefile and style file
switch p.Results.Dataset
    case 'plate'
        filename = '../../GIS/global_tectonics/plates&provinces/plates.shp';
        t = readtable('plate_type_styles.xlsx');
        field = {"plate","plate_type"};
        sheet = "";
    case 'province'
        filename = '../../GIS/global_tectonics/plates&provinces/global_gprv.shp';
        field = "prov_type";
        t = readtable('global_gprv_styles.xlsx','Sheet',field);
    case 'lastorogen'
        filename = '../../GIS/global_tectonics/plates&provinces/global_gprv.shp';
        field = "lastorogen";
        t = readtable('global_gprv_styles.xlsx','Sheet',field);
end
s = shaperead(filename);

t.color = [t.r t.g t.b]/255;

% setup or hold on axes
ax = gca;
hold(ax,'on');
prepmap;

cp = 1;
ca = 1;
for i = 1:length(s)
    if isempty(s(i).Y)
        continue;
    end

    x = [NaN; s(i).X'];
    y = [NaN; s(i).Y'];

    if length(field) == 1
        j = find(strcmp(t{:,field},getfield(s(i),field)));
    else
        j = find(strcmp(t{:,field{1}},getfield(s(i),field{1})) & strcmp(t{:,field{2}},getfield(s(i),field{2})));
    end
    
    ind = find(isnan(x));
    for k = 2:length(ind)
        for n = 1:4
            platepoly = polyshape([x(ind(k-1)+1:ind(k)-1); x(ind(k-1)+1)], [y(ind(k-1)+1:ind(k)-1); y(ind(k-1)+1)],'Simplify',true);
    
            % find intersection of polygons
            polyi = intersect(platepoly,quadpoly(n));
    
            if ~isempty(polyi.Vertices) & ~isempty(j)
                sp(cp) = s(i);

                % split lines at map octahedral boundaries
                [sp(cp).X,sp(cp).Y,xa,ya] = OctahedralProjection(polyi.Vertices(:,2),polyi.Vertices(:,1), ...
                    'Projection',p.Results.Projection, 'Style',p.Results.Style, ...
                    'Shift',p.Results.Shift, 'Inset',p.Results.Inset);

                polyi.Vertices = [sp(cp).X sp(cp).Y];

                plot(ax,polyi,'FaceColor',t.color(j,:),'FaceAlpha',0.4,'LineWidth',0.25,'EdgeColor','none');
                
                cp = cp + 1;
            end
        end
    end

    if p.Results.Inset && ~isempty(xa)
        xa = [NaN; xa];
        ya = [NaN; ya];

        ind = find(isnan(xa));
        for k = 2:length(ind)
            platepoly = polyshape([xa(ind(k-1)+1:ind(k)-1); xa(ind(k-1)+1)], [ya(ind(k-1)+1:ind(k)-1); ya(ind(k-1)+1)]);

            % find intersection of polygons
            polyi = intersect(platepoly,polyinset);
        end
        
        if ~isempty(polyi.Vertices) & ~isempty(j)
            plot(ax,polyi,'FaceColor',t.color(j,:),'FaceAlpha',0.4,'LineWidth',0.25,'EdgeColor','none');
    
            sa(ca) = s(i);
            sa(ca).X = polyi.Vertices(:,1);
            sa(ca).Y = polyi.Vertices(:,2);
    
            ca = ca + 1;
        end
    end
end

% write output shapefiles
filename = [filename,'_',p.Results.Projection,'_',p.Results.Style,'_Shift_',num2str(p.Results.Shift)];

if ~exist([filename,'.shp'],'file')
    shapewrite(sp,[filename,'.shp']);
    if p.Results.Inset
        shapewrite(sa,[filename,'_Inset.shp']);
    end
end

return