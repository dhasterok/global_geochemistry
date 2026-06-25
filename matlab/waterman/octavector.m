function [x,y,dx,dy,varargout] = octavector(lat,lon,de,dn,varargin)

% parse inputs
p = inputParser;
addRequired(p,'Latitude',@isnumeric);
addRequired(p,'Longitude',@isnumeric);
addRequired(p,'DeltaE',@isnumeric);
addRequired(p,'DeltaN',@isnumeric)
addParameter(p,'Projection','Waterman',@ischar);
addParameter(p,'Style','Butterfly',@ischar);
addParameter(p,'Shift',0,@isnumeric);
addParameter(p,'Transform','forward',@ischar);
addParameter(p,'Inset',false,@islogical);
addParameter(p,'MaxLat',-60,@isnumeric);
addParameter(p,'Save','',@ischar);

parse(p,lat,lon,de,dn,varargin{:});

[x,y] = OctahedralProjection(lat,lon,'Projection',p.Results.Projection, ...
    'Style',p.Results.Style, ...
    'Shift',p.Results.Shift);

[X1,Y1] = OctahedralProjection([lat lat lat],[lon-1e-4 lon lon+1e-4], ...
    'Projection',p.Results.Projection, ...
    'Style',p.Results.Style, ...
    'Shift',p.Results.Shift);

[X2,Y2] = OctahedralProjection([lat-1e-4 lat lat+1e-4],[lon lon lon], ...
    'Projection',p.Results.Projection, ...
    'Style',p.Results.Style, ...
    'Shift',p.Results.Shift);

[dx,dy] = octarotate(de,dn,X1,Y1,X2,Y2);

quiver(x,y,dx,dy);
hold on;

if p.Results.Inset
    ind = lat < p.Results.MaxLat;
    
    [~,~,xa,ya] = OctahedralProjection(lat(ind),lon(ind),'Projection',p.Results.Projection, ...
        'Style',p.Results.Style, ...
        'Shift',p.Results.Shift, ...
        'Inset',p.Results.Inset, ...
        'MaxLat',p.Results.MaxLat);
    
    [~,~,Xa1,Ya1] = OctahedralProjection([lat(ind) lat(ind) lat(ind)],[lon(ind)-1e-4 lon(ind) lon(ind)+1e-4], ...
        'Projection',p.Results.Projection, ...
        'Style',p.Results.Style, ...
        'Inset',p.Results.Inset, ...
        'MaxLat',p.Results.MaxLat);

    Xa1 = reshape(Xa1,length(lat),3);
    Ya1 = reshape(Ya1,length(lat),3);
    
    [~,~,Xa2,Ya2] = OctahedralProjection([lat(ind)-1e-4 lat(ind) lat(ind)+1e-4],[lon(ind) lon(ind) lon(ind)], ...
        'Projection',p.Results.Projection, ...
        'Style',p.Results.Style, ...
        'Inset',p.Results.Inset, ...
        'MaxLat',p.Results.MaxLat);

    Xa2 = reshape(Xa2,length(lat),3);
    Ya2 = reshape(Ya2,length(lat),3);
    
    [dxa,dya] = octarotate(de,dn,Xa1,Ya1,Xa2,Ya2);

    varargout{1} = xa;
    varargout{2} = ya;
    varargout{3} = dxa;
    varargout{4} = dya;

    quiver(xa,ya,dxa,dya);
end

prepmap

return


function [dx,dy] = octarotate(de,dn,X1,Y1,X2,Y2)

alpha = sum(atan2(diff(Y1,1,2),diff(X1,1,2)),2)/2;
beta = sum(atan2(diff(Y2,1,2),diff(X2,1,2)),2)/2;

dx = de.*cos(alpha) + dn.*cos(beta);
dy = de.*sin(alpha) + dn.*sin(beta);

return