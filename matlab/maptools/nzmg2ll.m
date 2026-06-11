function [lon,lat] = utm2ll(x,y);
% UTM2LL - Coverts UTM to Lat/Lon coordinates.
%
%    [LON,LAT] = UTM2LL(E,N,ZONE,DATUM) coverts easting, E, and
%    northing, N, from universal transverse mercator (UTM)
%    coordinates to longitude, LON, and latitude, LAT.   UTM2LL
%    requires the UTM zone, ZONE, and the DATUM for the reference
%    ellipsoid.
%
% Last Modified: 30 Sept. 2007

[a,b] = getdatum('nzmg');
ecc = eccentricity(a,b);

% Set false origin
%feast = 2510000;
%fnorth = 6023150;
feast = 1.6e6;
fnorth = 1e7;


b1 = 1.3231270439;
b2 = -0.577245789 - i*0.007809598;
b3 = 0.508307513 - i*0.112208952;
b4 = -.15094762 + i*0.18200602;
b5 = 1.01418179 + i*1.64497696;
b6 = 1.9660549 + i*2.5127645;

B1 = 0.7557853228;
B2 = 0.249204646 + i*0.003371507;
B3 = -0.001541739 + i*0.041058560;
B4 = -0.10162907 + i*0.01727609;
B5 = -0.22623489 - i*0.36249218;
B6 = -0.6870983 - i*1.1651967;

z = (y - fnorth + i*(x - feast))/a;

% First approximation of zeta
zeta = b1*z + b2*z.^2 + b3*z.^3 + b4*z.^4 + b5*z.^5 + b6*z.^6;

% Second approximation of zeta
zeta = (z + B2*zeta.^2 + 2*B3*zeta.^3 + 3*B4*zeta.^4 + 4*B5*zeta.^5 + 5*B6*zeta.^6)./ ...
       (B1 + 2*B2*zeta + 3*B3*zeta.^2 + 4*B4*zeta.^3 + 5*B5*zeta.^4 + 6*B6*zeta.^5);

dilat = real(zeta);

% Convert from isometric latitude to geographic latitude
% (note that it is in seconds and must be converted to degrees)

dlat = 1e5*(1.5627014243*dilat + 0.5185406398*dilat.^2 - 0.03333098*dilat.^3 ...
       - 0.1052906*dilat.^4 - 0.0368594*dilat.^5 + 0.007317*dilat.^6 ...
       + 0.01220*dilat.^7 + 0.00394*dilat.^8 - 0.0013*dilat.^9);

lat = -41 + dlat/3600;

% The computed longitude is in radians

lon = 173 + 180/pi*imag(zeta);

return
