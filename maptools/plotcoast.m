function c = plotcoast(varargin);
% PLOTCOAST - Plots coastline.
%
%    C = PLOTCOAST plots global coastlines and returns graphics
%    handles to coastline plot.
%
% 11 May 2011 By D. Hasterok (SIO)

C = [0 0 0];
lon0 = -180;

opt = 1;
if nargin > 0
    while opt < nargin
        switch lower(varargin{opt})
            case 'color'
                C = varargin{opt+1};
                opt = opt + 2;
            case 'lon0'
                lon0 = varargin{opt+1};
                opt = opt + 2;
            otherwise
                error('Unknown option.');
        end
    end
end

hold on;

load('coast.mat');
if lon0 == -180
    c = plot(clon,clat,'-','Color',C,'LineWidth',0.24);
else
    if lon0 < -180
        ind = clon > 360 + lon0;
        clon(ind) = clon(ind) - 360;
    else
        ind = clon < lon0;
        clon(ind) = clon(ind) + 360;
    end
    c = 1;
    for i = 1:length(ind)-1
        if (ind(i) == 1 & ind(i+1) == 0) || (ind(i) == 0 & ind(i+1) == 1)
            clon2([c:c+1],1) = [clon(i); NaN];
            clat2([c:c+1],1) = [clat(i); NaN];
            c = c + 2;
        else
            clon2(c,1) = clon(i);
            clat2(c,1) = clat(i);
            c = c + 1;
        end
    end
    c = plot(clon2,clat2,'-','Color',C,'LineWidth',0.24);
end


xlabel('Longitude');
ylabel('Latitude');

set(gca,'Box','on');
axis equal;
axis([lon0 lon0+360 -90 90]);

return
