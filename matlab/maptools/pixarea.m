function A = pixarea(lat1,lat2,lon1,lon2,varargin);
% PIXAREA - Area on a sphere.
%
%   A = pixarea(lat1,lat2,lon1,lon2) computes the fractional area of
%   quad point pixels on a sphere.
%
%   A = pixarea(lat1,lat2,lon1,lon2,geoid) computes the fractional area
%   using the reference ellipsoid, GEOID.

% Changelog:
%   Original: 3 Apr. 2012 by D. Hasterok

% convert from degrees to radians
lat1 = lat1*pi/180;
lon1 = lon1*pi/180;

lat2 = lat2*pi/180;
lon2 = lon2*pi/180;

% convert latitude to authalic latitude if a geoid is given
if nargin == 5
    geoid = varargin{1};
    lat1 = authalic(lat1,geoid);
    lat2 = authalic(lat2,geoid);
end

% compute area as fraction of sphere
A = 0.25 * abs(lon1 - lon2) .* abs(sin(lat1) - sin(lat2)) / pi;

return
