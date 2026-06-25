function [x,y,Z] = geo2octahedral_raster(filename,varargin)

addpath ../worldgrid
addpath ../maptools

% parse inputs
p = inputParser;

addParameter(p,'Projection','Waterman',@ischar);
addParameter(p,'Style','Butterfly',@ischar);
addParameter(p,'Shift',0,@isnumeric);
addParameter(p,'Inset',false,@islogical);
addParameter(p,'MaxLat',-60,@isnumeric);
addParameter(p,'Export',false,@islogical);

parse(p,varargin{:});

% load file and set up lat and lon
switch filename
    case {'../worldgrid/ETOPO1_Bed_g_gdal.grd','../worldgrid/ETOPO1_Ice_g_gdal.grd','../worldgrid/age.2020.1.GTS2012.1m.nc'}
        n = 180*60;

        lon = linspace(-180,180,2*n+1);
        lat = linspace(-90,90,n+1)';

        %z = ncread(filename,'z');
        z = reshape(ncread(filename,'z'),length(lon),length(lat))';

        if contains(filename,'ETOPO1')
            z = flipud(z);
        end
    case '../../GIS/global_gravity/global_grav.nc'
        n = 180*30;

        lon = linspace(-180,180,2*n+1);
        lat = linspace(-90,90,n+1)';

        z = reshape(ncread(filename,'z'),length(lon),length(lat))';
    case '../worldgrid/gpw-v4-population-density-adjusted-to-2015-unwpp-country-totals-rev11_2020_2pt5_min_asc/gpw_density_2020.nc'
        z = ncread(filename,'Band1')';

        n = min(size(z));
        dl = 180/n;
        lon = linspace(-180+dl/2,180-dl/2,2*n);
        lat = linspace(-90+dl/2,90-dl/2,n);
end

%lat = linspace(-90,90,180*30+1)';
%lon = linspace(-180,180,360*30+1);

%filename = '../worldgrid/MohoDepthSzwillus.nc';
%elev = ncread(filename,'crustalthickness');
%filename = '../worldgrid/age.2020.1.GTS2012.1m.nc';
%filename = ;
%elev = reshape(ncread(filename,'z'),length(lon),length(lat))';

%imagesc(lon,lat,z);
%return

[lon,lat] = meshgrid(lon,lat);

[x,y,Z] = octaraster(lat,lon,z,'Projection',p.Results.Projection, ...
    'Style',p.Results.Style, ...
    'Shift',p.Results.Shift, ...
    'Inset',p.Results.Inset, ...
    'MaxLat',p.Results.MaxLat);

figure;
imagesc(x,y,Z,'AlphaData',~isnan(Z));
hold on;
axis xy;
prepmap;
%octamask('Inset',p.Results.Inset,'MaxLat',p.Results.MaxLat);
octagraticule(15,'Shift',p.Results.Shift,'Color',[0.2 0.2 0.2],'Inset',p.Results.Inset,'MaxLat',p.Results.MaxLat);
octaborder('Inset',p.Results.Inset,'MaxLat',p.Results.MaxLat);

%wiki2 = load('wiki2_palette.csv')/255;
%centermap(wiki2,10,0,[-8000 6500]);

% if export requested
if p.Results.Export
    outfile = [filename(1:end-3),'_',p.Results.Projection,'_',p.Results.Style];
    xyz2ncdf(outfile,{'x','y','z'},{'none','none',''},'Data',x,y,Z');
end

return