function [alat,Re] = authalic(glat,geoid);
% AUTHALIC - Computes authalic latitude.
%
%   [alat,Re] = authalic(glat,geoid) converts the geodetic latitude,
%   GLAT, to authalic latitude, ALAT, using the reference geoid, GEOID.
%   
%   This transformation converts an ellipsoid into a sphere of
%   equivalent equal surface area.

[nr,nc] = size(glat);
glat = glat(:);

% Changelog:
%   Original: 3 Apr. 2012 by D. Hasterok

% Reference: Wolfram Mathworld, 2012,
% http://mathworld.wolfram.com/AuthalicLatitude.html

% Get ellipsoid parameters
[ecc,a,b] = eccentricity(geoid);

% convert geodetic latitude to radians
phi = glat*pi/180;

% compute scale factors
q = scalefactor(ecc,phi);
qp = scalefactor(ecc,0.5*pi);

% compute authalic latitude (in radians)
beta = asin(q/qp);

% convert authalic latidude to degrees
alat = 180/pi*beta;

alat = reshape(alat,nr,nc);

% Radius of sphere approximating ellipsoid
Re = a*sqrt(0.5*qp);

return


% scale factor for authalic latitude
function q = scalefactor(ecc,phi);

s = sin(phi);
q = (1 - ecc^2)*( s./(1 - ecc^2*s.^2) - 0.5/ecc*log((1 - ecc*s)./(1 + ecc*s)) );

return
