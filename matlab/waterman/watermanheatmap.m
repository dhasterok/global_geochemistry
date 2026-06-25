function watermanheatmap(data,field,varargin)

p = inputParser;
addParameter(p,'Projection','Waterman',@ischar);
addParameter(p,'Style','Butterfly',@ischar);
addParameter(p,'Inset',false,@islogical);
addParameter(p,'MaxLat',-60,@isnumeric);
addParameter(p,'Shift',0,@isnumeric);
addParameter(p,'Scale','linear',@ischar);

parse(p,varargin{:});


[X,Y] = OctahedralProjection(data.latitude,data.longitude,'Projection',p.Results.Projection, ...
    'Style',p.Results.Style, ...
    'Shift',p.Results.Shift);

switch p.Results.Scale
    case 'linear'
        Z = data.(field);
    case 'log'
        Z = log10(data.(field));
    otherwise
        error('Scale must be linear or log.')
end

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

% figure;
% plot(xc(~in),yc(~in),'b.');
% hold on;
% plot(xc(in),yc(in),'r.');
% axis equal

% create hexagon verticies (first is repeated to close polygon)
xv = dx*cos(rot) - dy*sin(rot);
yv = dx*sin(rot) + dy*cos(rot);

ind = find(~in);

%[x,y,Z] = geo2octahedral_raster('../worldgrid/ETOPO1_Ice_g_gdal.grd','Shift',20,'Inset',true);
%close all;

% create hexbin structure
% figure;
% hold on;
tempbin.xc = nan;
tempbin.yc = nan;
tempbin.xv = zeros(size(xv));
tempbin.yv = zeros(size(yv));
tempbin.n = 0;
tempbin.mean = nan;
tempbin.std = nan;

hexbin = repmat(tempbin,1,sum(~in));
for i = 1:sum(~in)
    hexbin(i).xc = xc(ind(i));
    hexbin(i).yc = yc(ind(i));
    hexbin(i).xv = xv + xc(ind(i));
    hexbin(i).yv = yv + yc(ind(i));

    % plot(hexbin(i).xv,hexbin(i).yv,'k-');
end
% axis equal

% compute the average of the dataset
h = waitbar(0,'Computing means...');
N = length(hexbin);
for i = 1:length(hexbin)
    if mod(i,600)
        waitbar(i/length(hexbin),h)
    end
    ind = min(hexbin(i).xv) <= X & X < max(hexbin(i).xv) & min(hexbin(i).yv) <= Y & Y < max(hexbin(i).yv);
    xtemp = X(ind);
    ytemp = Y(ind);
    ztemp = Z(ind);

    in = inpolygon(xtemp,ytemp,hexbin(i).xv,hexbin(i).yv);
    if sum(in) == 0
        continue;
    end
    n = sum(in);
    hexbin(i).n = n;
    hexbin(i).mean = mean(ztemp(in));
    if n > 1
        hexbin(i).std = std(ztemp(in));
    end
end
close(h)

figure;
hold on;
for i = 1:length(hexbin)
    % if hexbin(i).xc > -5.7 && hexbin(i).xc < 5.7 ...
    %     && ~(hexbin(i).xc < -2.3 && hexbin(i).yc < -2.75) && ~(hexbin(i).xc > 2.3 && hexbin(i).yc < -2.75)
    if ~isnan(hexbin(i).mean)
        fill(hexbin(i).xv,hexbin(i).yv,hexbin(i).mean,'LineWidth',0.1)
    end
end

octacoast('Shift',p.Results.Shift,'Inset',false);
%octaplatebndry('Shift',20,'Inset',false);
octaborder('Inset',false);

cb = colorbar;

colormap2('plasma');

pos = cb.Position;
cb.Position = [pos(1)-.1, pos(2)+.12, pos(3), pos(4)/3];

cb.Label.String = 'Heat Production (\muW m^{-3})';

cb.TickDirection = 'out';

clim([-1 log10(3)]);
cb.Ticks = log10([0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 2 3]);
cb.TickLabels = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 2 3];

return