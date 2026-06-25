addpath ../maptools

%im = imread('world.200408.3x5400x2700.png');
%dl = 360/5400;
%im = imread('world.200408.3x21600x10800.png');
%dl = 360/21600;
im = imread('../worldgrid/BlackMarble_2016_3km_geo.tif');
dl = 360/13500;

im = flipud(im);

imlat = linspace(-90+0.5*dl,90-0.5*dl,180/dl);
imlon = linspace(-180+0.5*dl,180-0.5*dl,360/dl);

[imlon,imlat] = meshgrid(imlon,imlat);

[x,y,Z] = octaraster(imlat,imlon,im,'Shift',20,'Inset',true,'MaxLat',-50);

imagesc(x,y,Z);
prepmap
axis xy
octaborder('Inset',true);