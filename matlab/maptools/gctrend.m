function varargout = gctrend(lon,lat,v,p1,p2,varargin)
% GCTREND - produces a 2D profile through a gridded dataset
%
%   [D,V] = gctrend(lon,lat,v,p1,p2) interpolates a grid to a profile line
%   from coordinates p1 to p2 line where D is the distance along the line 
%   and V is the interpolated value.
%
%   [plon,plat,V] = gctrend(lon,lat,v,p1,p2) returns the coordinate points
%   along a line from p1 to p2.
%
%   ... = gctrend(lon,lat,v,p1,p2,R) will change the discretization along
%   the interpolated line by a factor of R.
%
% 26 Feb 2022 by D. Hasterok

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

V = interp2(lon,lat,v,plon,plat);

switch nargout
    case 0
        return;
    case 2
        dx = C.Rearth * sphangle(plon,plat,p1(1),p1(2));
        
        varargout{1} = dx;
        varargout{2} = V;
    case 3
        varargout{1} = plon;
        varargout{2} = plat;
        varargout{3} = V;
    otherwise
        error('Incorrect number of arguments requested');
end

return