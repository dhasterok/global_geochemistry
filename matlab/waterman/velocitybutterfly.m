function velocitybutterfly(varargin)

addpath ../plotting
addpath ../toolbox
addpath ../maptools

p = inputParser;

addParameter(p,'Projection','Waterman',@ischar);
addParameter(p,'Style','Butterfly',@ischar);
addParameter(p,'Shift',0,@isnumeric);
addParameter(p,'Inset',false,@islogical);
addParameter(p,'MaxLat',-60,@isnumeric);
addParameter(p,'Export',false,@islogical);
addParameter(p,'InterpMethod','nearest',@ischar);
addParameter(p,'ResampleFactor',1,@isnumeric);

parse(p,varargin{:});
filename = '/Users/derrickhasterok/Google Drive/toolbox/geophysics/seismic/tomography/3D2018-08Sv-depth.nc';
lon = ncread(filename,'longitude');
lat = ncread(filename,'latitude');
depth = ncread(filename,'depth');

i = find(depth == 70);
z = ncread(filename,'dvs');

z = z(:,:,i)';

lon = [-179.99; lon];
lat = [lat; -90];

z = [0.5*(z(:,1) + z(:,end)) z];
z = [z; z(end,:)];

% filename = '../worldgrid/MohoDepthSzwillus.nc';
% lon = ncread(filename,'lon');
% lat = ncread(filename,'lat');
%  
% z = ncread(filename,'crustalthickness')';

[lon,lat] = meshgrid(lon,lat);

[x,y,Z] = octaraster(lat,lon,z, ...
    'Inset',p.Results.Inset, ...
    'Shift',p.Results.Shift, ...
    'MaxLat',p.Results.MaxLat, ...
    'InterpMethod',p.Results.InterpMethod, ...
    'ResampleFactor',p.Results.ResampleFactor);

[xb,yb,xbi,ybi] = panelborders(p);
[X,Y] = meshgrid(x,y);

ind = ~isnan(Z);
F = scatteredInterpolant(X(ind),Y(ind),Z(ind));

Z = F(X,Y);
% for j = 1:size(xb,2)
%     in = inpolygon(X,Y,xb(:,j),yb(:,j));
% 
%     if ~any(isnan(Z))
%         continue;
%     end
% 
%     ind = in & ~isnan(Z);
%     F = scatteredInterpolant(X(ind),Y(ind),Z(ind));
% 
%     ind = in & isnan(Z);
%     Z(ind) = F(X(ind),Y(ind));
% end
% 
% in = inpolygon(X,Y,xbi,ybi);
% ind = in & ~isnan(Z);
% F = scatteredInterpolant(X(ind),Y(ind),Z(ind));
% 
% ind = in & isnan(Z);
% Z(ind) = F(X(ind),Y(ind));

figure;
imagesc(x,y,Z);
prepmap
octacoast('Inset',p.Results.Inset,'Shift',p.Results.Shift,'MaxLat',p.Results.MaxLat);
octaborder('Inset',p.Results.Inset,'MaxLat',p.Results.MaxLat);
colormap2('rwb');

if p.Results.Export
    xyz2ncdf(['shear_velocity_',num2str(depth(i)),'km.nc'],{'x','y','dvs'},{'none','none','km/s'},'Data',x,y,Z');
    %xyz2ncdf('crustal_thickness',{'x','y','z'},{'none','none','km'},'Data',x,y,Z');
end

% write geotiff
geotiffwrite(['shear_velocity_',num2str(depth(i)),'km.tif'],Z, ...
    georasterref('RasterSize',size(Z), ...
    'LatitudeLimits',[min(y) max(y)], ...
    'LongitudeLimits',[min(x) max(x)]));

return