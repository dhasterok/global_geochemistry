function [A,lat] = gridwts(maxlat,dlat,dlon,varargin);
% GRIDWTS - Fractional area of square pixels on a sphere
%
%   [A,lat] = gridwts(maxlat,dlat,dlon) computes the fractional area of
%   a geographic pixelated grid with width DLON and height DLAT.  The
%   maximum latitude, MAXLAT, of the grid is required in order to
%   determine if the top pixel is at the pole.
%
%   [A,lat] = gridwts(maxlat,dlat,dlon,geoid) computes the fractional
%   areas on an ellipsoid, GEOID.

% Changelog:
%   Original: 3 Apr. 2012 by D. Hasterok

% Correction factor for grid cells that lie over the pole.
if maxlat == 90
    n = 2;
else
    n = 1;
end

% Pixel width
lon = [-dlon/2 dlon/2];

% Grid latitudes.  Since the corrections are symetric about the equator,
% only the northern hemisphere corrections are computed
lat = [maxlat-90:dlat:maxlat] - dlat/2;
lat = [lat 90];

if nargin == 4
    geoid = varargin{1};
end

% Pixel areas
A = pixarea(lat(1:end-1),lat(2:end),lon(1),lon(2),geoid);

% Correct area at pole
A(end) = n*A(end);

% Extend from -90 to 90
% Cell latitudes
lat = midpt(lat);
lat(end) = 90;
lat = [-fliplr(lat(2:end)) lat];

% Area
A = [fliplr(A(2:end)) A];

return
