function delta = sphangle(lon,lat,lon0,lat0,varargin);
% SPHANGLE - Spherical distance.
%
%    DELTA = SPHANGLE(LON,LAT,LON0,LAT0) computes the angle, DELTA,
%    between the set of points given by (LON,LAT) and a reference
%    (LON0,LAT0) on a sphere.  The angle DELTA is given in degrees.
%
%    DELTA = SPHANGLE(LON,LAT,LON0,LAT0,H) computes the signed angle,
%    DELTA, determined by the hemisphere, H.
%
%    Case H
%       0: (default) No hemisphere, all angles positive.
%       1: Hemisphere LON > LON0 positive
%       2: Hemisphere LON < LON0 positive
%       3: Hemisphere LAT > LAT0 positive
%       4: Hemisphere LAT < LAT0 positive
%
% Last Modified: 23 Sept. 2007 by D. Hasterok

% Get hemisphere if given, otherwise assume all angles will be positive
if nargin < 4 | nargin > 5
    error('ERROR (sphangle.m): Incorrect number of arguments.');
elseif nargin == 4
    h = 0;
else
    h = varargin{1};
end

% Convert to radians
theta0 = pi/180*lon0;
phi0 = pi/180*lat0;
% Compute reference in cartesian coordinates
[x0,y0,z0] = sph2cart(theta0,phi0,1);
p0 = [x0,y0,z0];

% Convert to radians
theta = pi/180*lon;
phi = pi/180*lat;
% Compute points to cartesian coordinats
r = ones([length(phi),1]);
[x,y,z] = sph2cart(theta,phi,r);
p = [x,y,z];

% Compute the angle between the vectors in degrees
ratio = dotprod(p0,p)./(normv(p0).*normv(p));
ratio(ratio > 1) = 1;
ratio(ratio < -1) = -1;

delta = acos(ratio);

% Set direction to delta if necessary
sgn = ones(size(theta));

switch h
    case 0
        return
    case 1
        change = find(theta < theta0);
    case 2
        change = find(theta > theta0);
    case 3
        change = find(phi < phi0);
    case 4
        change = find(phi > phi0);
    otherwise
        str = ['WARNING (sphangle.m): H must be integer from 0 to 4.';
               '   Assuming no hemisphere chosen.                   '];
        warning(str);
        return
end

sgn(change) = -1;
delta = delta.*sgn;

return


% Dot product of the vector with the reference
function d = dotprod(p0,p)

d = p0(1)*p(:,1) + p0(2)*p(:,2) + p0(3)*p(:,3);

return


function n = normv(v)
% Norm of the vector

n = sqrt(sum(v.^2,2));

return
