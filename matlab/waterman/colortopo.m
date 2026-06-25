function colortopo(varargin)

style = 'day';
if nargin == 1
    style = varargin{1};
end

switch style
    case 'day'
        % daytime mimicking satellite, "blue marble"
        brown = [0 177 143 87;
            1 245 244 242];
        green = [0 133 161 117;
            1 233 321 187];
        blue = [0 49 76 95;
            0.06 78 120 151;
            0.8 166 213 241;
            1 203 235 240];
        white = [0 147 185 226;
            1 250 250 250];
    case 'night'
        % nightime mimicking satellite, "black marble"
        brown = [0 50 47 76;
            1 43 40 64];
        green = [0 27 31 41;
            1 37 39 50];
        blue = [0 3 5 13;
            1 25 26 34];
        white = [0 42 40 84;
            1 63 60 111];
end

% colormap bounds
bounds.sea = [-8000 0];
bounds.ice = [0 4500];
bounds.land = [0 5000];
bounds.precip = [0 1000];
%bounds.precip = [0 20];

% load data
x = ncread('ETOPO1_Bed_g_gdal._Waterman_Butterfly.nc','X');
y = ncread('ETOPO1_Bed_g_gdal._Waterman_Butterfly.nc','Y');
elev_bed = ncread('ETOPO1_Bed_g_gdal._Waterman_Butterfly.nc','x')';

elev_ice = ncread('ETOPO1_Ice_g_gdal._Waterman_Butterfly.nc','x')';

xp = ncread('mean_annual_precipitation_Waterman_Butterfly_Shift_20_Inset.nc','X');
yp = ncread('mean_annual_precipitation_Waterman_Butterfly_Shift_20_Inset.nc','Y');

precip = ncread('mean_annual_precipitation_Waterman_Butterfly_Shift_20_Inset.nc','z')';

% interpolate precipitation to same positions as topography
[X,Y] = meshgrid(x,y);
[Xp,Yp] = meshgrid(xp,yp);

%imagesc(xp,yp,precip)
%precip = reshape(interp2(Xp,Yp,precip,x(:)',y(:)),length(y),length(x));
ind = ~isnan(precip);
F = scatteredInterpolant(Xp(ind),Yp(ind),precip(ind));
precip = reshape(F(X(:),Y(:)),length(y),length(x));

% normalize precipitation to bounds
precip = (precip - bounds.precip(1))/(bounds.precip(2) - bounds.precip(1));
precip(precip < 0) = 0;
precip(precip > 1) = 1;

% indicies for ice, sea, and land
ice = elev_ice - elev_bed > 0;
sea = elev_ice < 0;
land = ~ice & ~sea & ~isnan(elev_ice);

% create image matrix
im = repmat(255,length(y),length(x),3);

% compute color of sea
tmp = colorbyvalue(elev_ice(sea),bounds.sea,blue);
[i,j] = ind2sub(size(elev_ice),find(sea));

for k = 1:3
    ind = sub2ind([size(elev_ice) 3],i,j,repmat(k,length(i),1));
    im(ind) = tmp(:,k);
end

% compute color of ice
tmp = colorbyvalue(elev_ice(ice),bounds.ice,white);
[i,j] = ind2sub(size(elev_ice),find(ice));

for k = 1:3
    ind = sub2ind([size(elev_ice) 3],i,j,repmat(k,length(i),1));
    im(ind) = tmp(:,k);
end

% compute color of topography
% brown
btopo = colorbyvalue(elev_ice(land),bounds.land,brown);
% green
gtopo = colorbyvalue(elev_ice(land),bounds.land,green);

% mix by brown and green by relative precip (0 = brown, 1 = green)
tmp = repmat(precip(land),1,3).*gtopo + (1 - repmat(precip(land),1,3)).*btopo;
[i,j] = ind2sub(size(elev_ice),find(land));
for k = 1:3
    ind = sub2ind([size(elev_ice) 3],i,j,repmat(k,length(i),1));
    im(ind) = tmp(:,k);
end

% convert to uint8 for standard raster image
im = uint8(im);

figure;
imagesc(x,y,im);
prepmap
axis xy
octaborder('Inset',true,'MaxLat',-50);

% write geotiff
geotiffwrite(['global_topo_',style,'.tif'],im, ...
    georasterref('RasterSize',size(im), ...
    'LatitudeLimits',[min(y) max(y)], ...
    'LongitudeLimits',[min(x) max(x)]));

% legend
figure
n = 32;
leg = repmat(255,n,n,3);

subplot(1,4,2:3);
for i = 1:3
    leg(:,1,i) = linspace(brown(1,i+1),brown(2,i+1),n)';
    leg(:,n,i) = linspace(green(1,i+1),green(2,i+1),n)';

    for j = 1:n
        leg(j,:,i) = linspace(leg(j,1,i),leg(j,n,i),n);
    end
end
imagesc(linspace(bounds.precip(1),bounds.precip(2),n)/10,linspace(bounds.land(1),bounds.land(2),n)/1000,uint8(leg));
axis xy;
axis square;
title('Land');
xlabel('Precipitation (cm/a)');

subplot(144);
leg3 = repmat(255,n,1,3);
for i = 1:3
    leg3(:,1,i) = linspace(white(1,i+1),white(2,i+1),n);
end
imagesc(1,linspace(bounds.ice(2),bounds.ice(1),n)/1000,uint8(leg3));
title('Ice');
set(gca,'XTick',[],'TickDir','out');
axis xy;

n = 101;
subplot(141);
leg2 = repmat(255,n,1,3);
for i = 1:3
    for k = 2:size(blue,1)
        leg2(blue(k-1,1)*100+1:blue(k,1)*100+1,1,i) = linspace(blue(k-1,i+1),blue(k,i+1),blue(k,1)*100-blue(k-1,1)*100+1)';
    end
end
imagesc(1,linspace(bounds.sea(1),bounds.sea(2),n)/1000,uint8(leg2));
title('Water');
ylabel('Elevation (km)');
set(gca,'XTick',[],'TickDir','out');
axis xy;

return