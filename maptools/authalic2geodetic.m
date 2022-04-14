function [glat,Re] = authalic2geodetic(alat,geoid);

[nr,nc] = size(alat);
alat = alat(:);

% Get ellipsoid parameters
[ecc,a,b] = eccentricity(geoid);

% convert geodetic latitude to radians
phi = alat*pi/180;



% compute scale factors
q = scalefactor(ecc,phi);
qp = scalefactor(ecc,0.5*pi);

% compute authalic latitude (in radians)
beta = asin(q/qp);




% convert authalic latidude to degrees
glat = 180/pi*beta;

glat = reshape(glat,nr,nc);

% Radius of sphere approximating ellipsoid
Re = a*sqrt(0.5*qp);

return


% scale factor for authalic latitude
function q = scalefactor(ecc,phi);

s = sin(phi);
q = (1 - ecc^2)*( s./(1 - ecc^2*s.^2) - 0.5/ecc*log((1 - ecc*s)./(1 + ecc*s)) );

return