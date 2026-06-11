function [lon,lat] = utm2ll(x,y,Z,datum);
% UTM2LL - Coverts UTM to Lat/Lon coordinates.
%
%    [LON,LAT] = UTM2LL(E,N,ZONE,DATUM) coverts easting, E, and
%    northing, N, from universal transverse mercator (UTM)
%    coordinates to longitude, LON, and latitude, LAT.   UTM2LL
%    requires the UTM zone, ZONE, and the DATUM for the reference
%    ellipsoid.
%
%    If zone is negative, then the Northing is referenced to the south
%    pole, not the equator.
%
% Last Modified: 2 Feb. 2022
% fixed southern hemisphere latitude
datum = lower(datum);

% Set false origin
switch datum
    case {'australian1984'; 'gda84'}
        feast = 0;
        fnorth = 1e7;
    case {'australian1994', 'gda94'}
        feast = 5e5;
        fnorth = 1e7;
    case 'nzmg'
        feast = 2510000;
        fnorth = 6023150;
        txt = []';
        error('The NZMG is an orthomorphic projection and this script cannot be used. Use NZMG2LL instead.');
    case 'nztm'
        feast = 1600000;
        fnorth = 10000000;
    case 'ykj' % finland
        datum = 'intl1924';
        Z = 35;
        feast = 3500000;
        fnorth = 0;
    case 'rt90' % sweden
        fnorth = 0;
        feast = 1.5e6;
        lon0 = 15 + 48/60 + 29.8/3600;
        
        [a,b] = getdatum('bessel1841');
    otherwise
        feast = 5e5;
        fnorth = 0;
end

sgn = sign(Z);
Z = abs(Z);

if length(Z) == 1 & length(x) > 1
    sgn = sgn*ones(size(x));
end

    

if ~strcmp(datum,'rt90')
    [a,b] = getdatum(datum);
    lon0 = 6*pi/180*(Z - 1) + 3*pi/180 - pi;
end
ecc = eccentricity(a,b);

k0 = 0.9996;

x = x - feast;
y = y - fnorth;


e1 = (1 - (1 - ecc^2)^0.5)/(1 + (1 - ecc^2)^0.5);

%M0 = a*((1 ...
%           - ecc^2/4 - 3*ecc^4/64 - 5*ecc^6/256)*lat0 ...
%           - (3*ecc^2/8 + 3*ecc^4/32 + 45*ecc^6/1024)*sin(2*lat0) ...
%           + (15*ecc^4/256 + 45*ecc^6/1024)*sin(4*lat0) ...
%           - (35*ecc^6/3072)*sin(6*lat0));

       
%M = M0 + y/k0;
M = y/k0;

mu = M/(a*(1 - ecc^2/4 - 3*ecc^4/64 - 5*ecc^6/256 - 7*ecc^8/1024));

lat1 = mu + (3*e1/2 - 27*e1^3/32)*sin(2*mu) ...
          + (21*e1^2/16 - 55*e1^4/32)*sin(4*mu) ...
          + (151*e1^3/96)*sin(6*mu) ...
          + (1097*e1^4/512)*sin(8*mu);

ep = ecc^2/(1 - ecc^2);
C1 = ep^2*(cos(lat1)).^2;
T1 = (tan(lat1)).^2;
N1 = a./(1 - ecc^2*(sin(lat1)).^2).^0.5;
R1 = a*(1 - ecc^2)./(1 - ecc^2*(sin(lat1)).^2).^1.5;
D = x./(N1*k0);

lat = 180/pi*(lat1 ...
          - (N1.*tan(lat1)./R1)...
          .*(D.^2/2 - (5 + 3*T1 + 10*C1 - 4*C1.^2 - 9*ep^2).*D.^4/24 ...
          + (61 + 90*T1 + 298*C1 + 45*T1.^2 - 252*ep^2 - 3*C1.^2).*D.^6/720));

lat(sgn >= 0) = 90 - lat(sgn >= 0);
lat(sgn < 0) = -90 + lat(sgn < 0);

lon = 180/pi*(lon0 ...
          + (D - (1 + 2*T1 + C1).*D.^3/6 ...
          + (5 - 2*C1 + 28*T1 - 3*C1.^2 + 8*ep^2 + 24*T1.^2).*D.^5/120)./cos(lat1));

return
