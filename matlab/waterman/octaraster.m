function varargout = octaraster(glat,glon,z,varargin)
% OCTARASTER - converts a geographic raster to octahedral
%
%   h = octagraticule(lat,lon,z) converts a geographic raster to an
%   octahedral raster.  lat and lon should be produced by meshgrid
%
%   There are several option value pairs:
%
%       'Projection'        choose between Waterman and Cahill
%                           Cahill not yet implemented
%
%       'Style'             'Butterfly' or 'M'
%
%       'Shift'             shifts the map latitude to the east by the 
%                           specified angle in degrees
%
%       'Inset'             plot Antarctica as contiguous, true or false
%                           default is false
%
%       'MaxLat'            if an inset is included, then a maximum
%                           latitude extent for the inset may be included.
%                           default is -60, but -50 will capture the entire
%                           Scotia Plate
%
%       'ResampleFactor'    integer resampling factor
%
%       'InterpMethod'      interpolation method for resampling, default is
%                           'nearest'
%
% see also OCTAHEDRALPROJECTION, OCTABORDER, OCTACOAST, OCTARASTER


%parse inputs
p = inputParser;

addRequired(p,'Latitude',@isnumeric);
addRequired(p,'Longitude',@isnumeric);
addRequired(p,'Values',@isnumeric);
addParameter(p,'Projection','Waterman',@ischar);
addParameter(p,'Style','Butterfly',@ischar);
addParameter(p,'Shift',0,@isnumeric);
addParameter(p,'Inset',false,@islogical);
addParameter(p,'MaxLat',-60,@isnumeric);
addParameter(p,'ResampleFactor',1,@isnumeric);
addParameter(p,'InterpMethod','nearest',@ischar);

parse(p,glat,glon,z,varargin{:});

% determine resolution
dl = max(abs([glon(1,1)-glon(1,2), glon(1,1)-glon(2,1)]));

% comute border vertices
if p.Results.Inset
    [xb,yb,xa,ya] = panelborders(p);
else
    [xb,yb] = panelborders(p);
end

switch p.Results.Style
    case 'Butterfly'
        % Waterman Butterfly dimensions
        if p.Results.Inset
            % w/ Antarctica
            rng = [-6.5622    6.5622   -5.4898    3.8561];
        else
            % w/out Antarctica
            rng = [-6.5622    6.5622   -3.6830    3.8561];
        end
    case 'M'
        % Waterman M dimensions
        if p.Results.Inset
            % w/ Antarctica
            rng = [-6.9282    6.9282   -4.4898    2.6830];
        else
            % w/out Antarctica
            rng = [-6.9282    6.9282   -2.6830    2.6830];
        end
    otherwise
        error('Unknown style.');
end
nx = p.Results.ResampleFactor*360/dl;
ny = p.Results.ResampleFactor*round(nx*(rng(4) - rng(3))/(rng(2) - rng(1)));

% produce grid over map area
x = linspace(rng(1),rng(2),nx)';
y = linspace(rng(3),rng(4),ny);

% produce empty value matrix
[X,Y] = meshgrid(x,y);
image_flag = false;
if size(z,3) == 1
    Z = nan(size(X));
elseif size(z,3) == 3
    image_flag = true;
    z = double(z);
    Z = repmat(255,size(X,1),size(X,2),3);
else
    error('Raster must have on z value or be an RGB color triplet.');
end

% find only values within map area
ind = false(size(X));
ind = ind(:);
parfor i = 1:4
    [in,on] = inpolygon(X(:),Y(:),xb(:,i),yb(:,i));
    ind = ind | in | on;
end

% convert to lat,lon
[lat,lon] = OctahedralProjection(X(ind),Y(ind), ...
    'Projection',p.Results.Projection, 'Style',p.Results.Style, ...
    'Shift',p.Results.Shift, 'Transform','inverse');

% determine value
if image_flag
    [i,j] = ind2sub(size(X),find(ind));
    for k = 1:3
        ind = sub2ind(size(Z),i,j,repmat(k,size(i)));
        Z(ind) = interp2(glon,glat,z(:,:,k),lon,lat);
    end
else
    Z(ind) = interp2(glon,glat,z,lon,lat);
end

if p.Results.Inset
    [in,on] = inpolygon(X(:),Y(:),xa,ya);
    
    inda = in | on;

    [~,~,alat,alon] = OctahedralProjection(X(inda),Y(inda), ...
        'Projection',p.Results.Projection, ...
        'Style',p.Results.Style, ...
        'Transform','inverse', ...
        'Inset',true, ...
        'Shift',p.Results.Shift, ...
        'MaxLat', p.Results.MaxLat);
    
    [i,j] = ind2sub(size(X),find(inda));
    if image_flag
        for k = 1:3
            inda = sub2ind(size(Z),i,j,repmat(k,size(i)));
            Z(inda) = interp2(glon,glat,z(:,:,k),alon,alat,p.Results.InterpMethod);
        end
    else
        Z(inda) = interp2(glon,glat,z,alon,alat);
    end
end

if nargout == 1
   figure;
   ax = gca;
   hold(ax,'on');
   imagesc(ax,x,y,Z);

   varargout = ax;
elseif nargout == 3
   if image_flag
       Z = uint8(Z);
   end
   varargout = {x,y,Z};
end

return