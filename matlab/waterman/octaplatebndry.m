function varargout = octaplatebndry(varargin)

% parse inputs
p = inputParser;
addParameter(p,'Projection','Waterman',@ischar);
addParameter(p,'Style','Butterfly',@ischar);
addParameter(p,'Shift',0,@isnumeric);
addParameter(p,'Transform','forward',@ischar);
addParameter(p,'Inset',false,@islogical);
addParameter(p,'MaxLat',-50,@isnumeric);
addParameter(p,'Export',false,@islogical);
addParameter(p,'Scale',1,@isnumeric)

parse(p,varargin{:});

scale = p.Results.Scale;

% plate boundary styles
t = readtable('plate_boundary_styles.xlsx');
t.color = [t.r t.g t.b]/255;

% read plate boundaries
%s = shaperead('plates_Waterman_Butterfly_Shift_20_Inset.shp');
s = shaperead('../../GIS/global_tectonics/plates&provinces/boundaries.shp');

% initialize output shapefile
sb = s;

ax = gca;
hold(ax,'on');
prepmap;

c = 1;
for i = 1:length(s)
    if isempty(s(i).Y)
        continue;
    end
    % split lines at map octahedral boundaries
    [x,y,xa,ya] = cutpolyconvert(s(i).Y,s(i).X,p);
    x = scale*x;
    y = scale*y;
    xa = scale*xa;
    ya = scale*ya;

    % output shapefile
    sb(i).X = x';
    sb(i).Y = y';

    % inset shapefile
    if p.Results.Inset && ~isempty(xa)
        si(c) = s(i);
        si(c).X = xa';
        si(c).Y = ya';

        c = c + 1;
    end
    
    j = find(strcmp(s(i).type,t.type));
    if s(i).level == 1
        lw = 0.75;
    else
        lw = 0.375;
    end

    % plot boundary line
    switch t.style{j}
        case 'solid'
            plot(ax,x,y,'-','Color',t.color(j,:));
            
            if p.Results.Inset && ~isempty(xa)
                plot(ax,xa,ya,'-','Color',t.color(j,:),'LineWidth',lw);
            end
        case 'dash'
            plot(ax,x,y,'--','Color',t.color(j,:),'LineWidth',lw);
            if p.Results.Inset && ~isempty(xa)
                plot(ax,xa,ya,'--','Color',t.color(j,:),'LineWidth',lw);
            end
        case 'barb'
            barbline(x,y,0.1,'Color',t.color(j,:),'Scale',0.15,'LineWidth',lw);
            if p.Results.Inset && ~isempty(xa)
                barbline(xa,ya,0.1,'Color',t.color(j,:),'Scale',0.15,'LineWidth',lw)
            end
    end
end

% write output shapefiles
if p.Results.Export
    filename = ['boundaries_',p.Results.Projection,'_',p.Results.Style,'_Shift_',num2str(p.Results.Shift)];
    shapewrite(sb,[filename,'.shp']);

    if p.Results.Inset
        shapewrite(si,[filename,'_Inset.shp']);
    end
end

return