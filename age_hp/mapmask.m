worldmap world
load coast
plotm(clat,clon)

[lon,lat] = meshgrid(linspace(-179, 179, 180),...
    linspace(-89.5,89.5, 180));

test = data(data.latitude>=-90 & data.latitude<=90 & data.longitude >=-180 & data.longitude <=180 & ~isnan(data.latitude) & ~isnan(data.longitude) & ~isnan(data.heat_production),:);
test = test(1:10000,:);
h=plotm(test.latitude(1:10000,:),test.longitude(1:10000,:),'.r');

%Landmask
ind = ~landmask(test.latitude,test.longitude,100);

test.latitude(ind,:) = NaN;
test.longitude(ind,:) = NaN;
h=plotm(test.latitude(1:10000,:),test.longitude(1:10000,:),'.b');