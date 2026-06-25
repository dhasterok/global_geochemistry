addpath ../maptools

data = load('~/Downloads/cru_ts4.08.1901.1910.pre.dat');

lat = data(:,1);
lon = data(:,2);

dl = 1/2;
Lon = [-180+dl/2:dl:180-dl/2];
Lat = [-90+dl/2:dl:90-dl/2]';

precip = nan(length(Lat),length(Lon),12);
i = floor((lat + 90)/dl + 1);
j = floor((lon + 180)/dl + 1);

for k = 1:12
    s = sub2ind([length(Lat),length(Lon),12],i,j,repmat(k,length(lat),1));

    precip(s) = data(:,k+2);
end

annual_precip = sum(precip(:,:,:),3);

[LON,LAT] = meshgrid(Lon,Lat);
whos
[x,y,Z] = octaraster(LAT,LON,annual_precip,'Shift',20,'Inset',true,'MaxLat',-50);

figure;
imagesc(x,y,Z);
prepmap;
octaborder('Inset',true,'MaxLat',-50);

xyz2ncdf('mean_annual_precipitation_Waterman_Butterfly_Shift_20_Inset',{'x','y','z'},{'none','none',''},'Data',x,y,Z');