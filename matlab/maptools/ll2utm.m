function [x,y,zone] = ll2utm(lon,lat,datum);
% LL2UTM - Converts from Lat/Lon to UTM coordinates.
%
%    [E,N,ZONE] = LL2UTM(LON,LAT,DATUM) converts latitude, LAT, longitude,
%    LON, to universal transverse mercator (UTM) coordinates, easting, E,
%    northing, N, within zone, ZONE.  To convert to UTM, it requires the 
%    DATUM for the reference ellipsoid.
%
% Last Modified: 23 Apr. 2008

[a,b] = getdatum(datum);
ecc = eccentricity(a,b);

switch datum
    case 'australian1984'
        feast = 0;
        fnorth = 1e7;
    case 'australian1994'
        feast = 5e5;
        fnorth = 1e7;
    case 'nzmg'
        feast = 2510000;
        fnorth = 6023150;
        error('The NZMG is an orthomorphic projection and this script cannot be used. Use NZMG2LL instead.');
    case 'nztm'
        feast = 1600000;
        fnorth = 10000000;
    otherwise
        feast = 5e5;
        fnorth = 0;
end

ind = find(lon >= 180);
lon(ind) = lon(ind) - 360;

% Compute the UTM zone
zone = div(lon+180,6) + 1;
lon0 = 6*pi/180*(zone - 1) + 3*pi/180 - pi;

lat = pi/180*lat;
lon = pi/180*lon;

k0 = 0.9996;

ep = ecc^2/(1 - ecc^2);
N = a./sqrt(1 - ecc^2*(sin(lat)).^2);
T = (tan(lat)).^2;
C = ep*(cos(lat)).^2;
A = (lon - lon0).*cos(lat);
M = a*((1 - ecc^2/4 - 3*ecc^4/64 - 5*ecc^6/256)*lat ...
        - (3*ecc^2/8 + 3*ecc^4/32 + 45*ecc^6/1024)*sin(2*lat) ...
        + (15*ecc^4/256 + 45*ecc^6/1024)*sin(4*lat) ...
        - 35*ecc^6/3072*sin(6*lat));

% Easting
x = k0*N.*(A + (1 - T + C).*A.^3/6 ...
        + (5 - 18*T + T.^2 + 72*C - 58*ep).*A.^5/120);

x = x + feast;


% Northing
y = k0*(M + N.*tan(lat).*(A.^2/2 + (5 - T + 9*C + 4*C.^2) ...
        .* A.^4/24 + (61 - 58*T + T.^2 + 600*C - 300*ep).*A.^6/720));

y = y + fnorth;

return
