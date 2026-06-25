function varargout = octagraticule(dl,varargin)
% OCTAGRATICULE - produces graticules on an octahedral plot
%
%   h = octagraticule(dl) produces graticules on an octahedral map with
%   spacing dlº.  Graphics handles for the graticules are returned in the
%   array h.
%
%   There are several option value pairs:
%
%       'Projection'        choose between Waterman and Cahill
%                           Cahill not yet implemented
%
%       'Style'             'Butterfly' or 'M'
%
%       'Shift'             shifts the map latitude to the east by the 
%                           specified angle in degrees
%
%       'Inset'             plot Antarctica as contiguous, true or false
%                           default is false
%
%       'MaxLat'            if an inset is included, then a maximum
%                           latitude extent for the inset may be included.
%                           default is -60, but -50 will capture the entire
%                           Scotia Plate
%
%       'Export'            if true, will export as a shapefile
%
% see also OCTAHEDRALPROJECTION, OCTABORDER, OCTACOAST, OCTARASTER

% parse inputs
p = inputParser;
addRequired(p,'Step',@isnumeric);
addParameter(p,'Projection','Waterman',@ischar);
addParameter(p,'Style','Butterfly',@ischar);
addParameter(p,'Shift',0,@isnumeric);
addParameter(p,'Inset',false,@islogical);
addParameter(p,'Color',[0.7 0.7 0.7],@isnumeric);
addParameter(p,'MaxLat',-60,@isnumeric);
addParameter(p,'Export',false,@islogical);

parse(p,dl,varargin{:});

ax = gca;
hold(ax,'on');
prepmap;

% meridians
% need to shift
lon = [-180:dl:180];

% need to split at equator to prevent artifacts
c = 1;
lat = [-90+dl:-1 -0.0001 NaN 0:90-dl];
for i = 1:length(lon)
    if p.Results.Inset
        [xg,yg,xa,ya] = OctahedralProjection(lat,repmat(lon(i),1,length(lat)), ...
            'Projection',p.Results.Projection, ...
            'Style',p.Results.Style, ...
            'Shift',p.Results.Shift, ...
            'Inset',p.Results.Inset, ...
            'MaxLat',p.Results.MaxLat);
        
        s(c) = geoshape(yg,xg,'Geometry','line','Meridian',lon(i),'Parallel',NaN);
        s(c+1) = geoshape(ya,xa,'Geometry','line','Meridian',lon(i),'Parallel',NaN);
        c = c + 2;
        
        plot(ax,xa,ya,'Color',p.Results.Color,'LineWidth',0.1);
    else
        %  convert to octahedral coordinates
        [xg,yg] = OctahedralProjection(lat,repmat(lon(i),1,length(lat)), ...
            'Projection',p.Results.Projection, ...
            'Style',p.Results.Style, ...
            'Shift',p.Results.Shift);

        s(c) = geoshape(yg,xg,'Geometry','line','Meridian',lon(i),'Parallel',NaN);
        c = c + 1;
    end

    % plot graticule
    h(i) = plot(ax,xg,yg,'Color',p.Results.Color,'LineWidth',0.1);
end

% parallels
% don't need to shift
lon = [linspace(-179.9999,-90.0001,90);
    linspace(-89.9999,-0.0001,90);
    linspace(0.0001,89.9999,90);
    linspace(90.0001,179.9999,90)];
lat = [-90+dl:dl:0-dl -1e-4 NaN 0:dl:90-dl];

k = length(h);
for j = 1:4
    for i = 1:length(lat)
        if p.Results.Inset && lat(i) <= -60
            [xg,yg,xa,ya] = OctahedralProjection(repmat(lat(i),1,size(lon,2)),lon(j,:), ...
                'Projection',p.Results.Projection, ...
                'Style',p.Results.Style, ...
                'Inset',true, ...
                'MaxLat',p.Results.MaxLat);
            
            s(c) = geoshape(yg,xg,'Geometry','line','Meridian',NaN,'Parallel',lat(i));
            s(c+1) = geoshape(ya,xa,'Geometry','line','Meridian',NaN,'Parallel',lat(i));
            c = c + 2;

            % plot border
            plot(ax,xa,ya,'Color',p.Results.Color,'LineWidth',0.1);
        else
            % convert to octahedral coordinates
            [xg,yg] = OctahedralProjection(repmat(lat(i),1,size(lon,2)),lon(j,:), ...
                'Projection',p.Results.Projection, ...
                'Style',p.Results.Style);

            s(c) = geoshape(yg,xg,'Geometry','line','Meridian',NaN,'Parallel',lat(i),'Parallel',lat(i));
            c = c + 1;
        end

        % plot graticule
        h(k+(j-1)*length(lat)+i) = plot(ax,xg,yg,'Color',p.Results.Color,'LineWidth',0.1);
    end
end

if nargout == 1
    varargout{1} = h;
end

% if Export, save as shapefile
if p.Results.Export
    filename = ['graticules_',p.Results.Projection,'_',p.Results.Style,'_Step_',num2str(p.Results.Step),'_Shift_',num2str(p.Results.Shift),'.shp'];
    shapewrite(s,filename);
end

return