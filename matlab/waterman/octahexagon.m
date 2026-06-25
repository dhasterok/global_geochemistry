close all;
clear all;

tmp = load('waterman_butterfly_mask_inset.csv');

xv = tmp(:,1);
yv = tmp(:,2);

xv = [xv; -flipud(xv)];
yv = [yv; flipud(yv)];

% rotation angles
rot = linspace(0,360,7)*pi/180;

% create first vertex
%dy = 0.5/(3*n*tan(pi/6));
yc = linspace(-5.125,3.80,107);
%yc = linspace(-5.125,3.80,143);
dy = (yc(2) - yc(1))/2;
dx = dy*tan(pi/6);

Nx = ceil(6.5622/(6*dx));
xc = [-Nx:Nx]*6*dx;

[xc,yc] = meshgrid(xc,yc);

xc2 = xc + 3*dx;
yc2 = yc + dy;

xc = [xc(:); xc2(:)];
yc = [yc(:); yc2(:)];

in = inpolygon(xc,yc,xv,yv);

figure;
plot(xc(~in),yc(~in),'b.');
hold on;
plot(xc(in),yc(in),'r.');
axis equal

% create hexagon verticies (first is repeated to close polygon)
xv = dx*cos(rot) - dy*sin(rot);
yv = dx*sin(rot) + dy*cos(rot);


ind = find(~in);

%[x,y,Z] = geo2octahedral_raster('../worldgrid/ETOPO1_Ice_g_gdal.grd','Shift',20,'Inset',true);
%close all;


% create hexbin structure
figure;
hold on;
tempbin.xc = nan;
tempbin.yc = nan;
tempbin.xv = zeros(size(xv));
tempbin.yv = zeros(size(yv));
tempbin.mean = nan;
tempbin.land = nan;
tempbin.sea = nan;
hexbin = repmat(tempbin,1,sum(~in));
parfor i = 1:sum(~in)
    hexbin(i).xc = xc(ind(i));
    hexbin(i).yc = yc(ind(i));
    hexbin(i).xv = xv + xc(ind(i));
    hexbin(i).yv = yv + yc(ind(i));

    plot(hexbin(i).xv,hexbin(i).yv,'k-');
end
axis equal

%[X,Y] = meshgrid(x,y);

%ind = ~isnan(Z);
%X = X(ind);
%Y = Y(ind);
%Z2 = Z(ind);
%save('topo.mat','X','Y','Z2')
load('topo.mat');

h = waitbar(0,'Computing means...');
N = length(hexbin);
for i = 1:length(hexbin)
    if mod(i,600)
        waitbar(i/length(hexbin),h)
    end
    ind = min(hexbin(i).xv) <= X & X < max(hexbin(i).xv) & min(hexbin(i).yv) <= Y & Y < max(hexbin(i).yv);
    xtemp = X(ind);
    ytemp = Y(ind);
    ztemp = Z2(ind);

    in = inpolygon(xtemp,ytemp,hexbin(i).xv,hexbin(i).yv);

    hexbin(i).mean = mean(ztemp(in));
    hexbin(i).land = sum(ztemp(in)>=0);
    hexbin(i).sea = sum(ztemp(in)<0);
    %hexbin(i).gmean = geomean(Z2(in));
end
close(h)

figure;
hold on;
for i = 1:length(hexbin)
    if (hexbin(i).mean > 0 || hexbin(i).land/(hexbin(i).land + hexbin(i).sea) >= 0.5)
    %if hexbin(i).land/(hexbin(i).land + hexbin(i).sea) > 0.4
        fill(hexbin(i).xv,hexbin(i).yv,hexbin(i).mean)
    end
end

octacoast('Shift',20,'Inset',true);

scale = 152.9452;
figure;
hold on;
dh = 10/scale;
for i = 1:length(hexbin)
    if (hexbin(i).mean > 0 || hexbin(i).land/(hexbin(i).land + hexbin(i).sea) >= 0.5)
    %if hexbin(i).land/(hexbin(i).land + hexbin(i).sea) > 0.4
        h = hexbin(i).mean/1e4;
        if h < 0
            h = 0;
        end
        h = 0.5*h + dh;
        fillhex(hexbin(i).xv*scale, hexbin(i).yv*scale, h*scale)
    end
end
colormap(copper);
clim([0 50]);
axis equal;
view(3);
lighting gouraud;
material dull;
light('Position',[1 1 1], 'Style','infinite');

ind = ([hexbin.mean] > 0 | [hexbin.land]./([hexbin.land] + [hexbin.sea]) >= 0.5);

hexbin_scaled = hexbin(ind);
for i = 1:length(hexbin_scaled)
    hexbin_scaled(i).xc = hexbin_scaled(i).xc*scale;
    hexbin_scaled(i).yc = hexbin_scaled(i).yc*scale;
    hexbin_scaled(i).xv = hexbin_scaled(i).xv*scale;
    hexbin_scaled(i).yv = hexbin_scaled(i).yv*scale;
    
    h = hexbin_scaled(i).mean/1e4;
    if h < 0
        h = 0;
    end
    hexbin_scaled(i).mean = scale*(0.5*h + dh);
end

figure;
hold on;
for i = 1:length(hexbin_scaled)
    fill(hexbin_scaled(i).xv, hexbin_scaled(i).yv, [0.75 0.75 0.75])
    text(hexbin_scaled(i).xc, hexbin_scaled(i).yc, ...
        num2str(round(hexbin_scaled(i).mean,1)), ...
        'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',6);
end
axis equal;

octacoast('Shift',20,'Inset',true,'Scale',scale);
%octaplatebndry('Shift',20,'Inset',true,'Scale',scale);
octaborder('Inset',true,'Scale',scale);
