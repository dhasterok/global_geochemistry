function varargout = gcslice(lon,lat,z,v,p1,p2,varargin)
% GCSLICE - produces a 2D profile through a gridded dataset
%
%   [D,Z,V] = gcslice(lon,lat,v,p1,p2) interpolates a 3D grid to a 2D
%   cross-section from coordinates p1 to p2 line where D is the distance 
%   along the line and V is the interpolated value at depths Z.
%
%   [plon,plat,Z,V] = gcslice(lon,lat,v,p1,p2) returns the map coordinate 
%   points along a line from p1 to p2.
%
%   ... = gcslice(lon,lat,v,p1,p2,R) will change the discretization along
%   the interpolated line by a factor of R.
%
% 28 Feb 2022 by D. Hasterok

p = inputParser;
addOptional(p,'ResampleFactor',1,@isnumeric);

parse(p,varargin{:});

r = p.Results.ResampleFactor;

% Earth's Radius
C = Constants;

dl = abs(lon(2) - lon(1));

% distance between end points
d = 180/pi*sphangle(p1(1),p1(2),p2(1),p2(2));

[plon,plat] = gcpoints(p1,p2,ceil(r*d/dl));

plot(plon,plat,'.')

% resample depth
if r ~= 1
    dz = diff(z)/r;

    dz = repmat(dz',r,1);
    dz = dz(:);
    Z = z(1) + [0; cumsum(dz)];
else
    Z = z;
end
X = repmat(plon,length(Z),1);
Y = repmat(plat,length(Z),1);
ZZ = repmat(Z',length(plon),1);

V = reshape(interp3(lon,lat,z,v,X,Y,flipud(ZZ(:))),length(plon),length(Z))';

switch nargout
    case 0
        return;
    case 3
        dx = C.Rearth * sphangle(plon,plat,p1(1),p1(2));
        
        varargout{1} = dx;
        varargout{2} = Z;
        varargout{3} = V;
    case 4
        varargout{1} = plon;
        varargout{2} = plat;
        varargout{3} = Z;
        varargout{4} = V;
    otherwise
        error('Incorrect number of arguments requested');
end

return