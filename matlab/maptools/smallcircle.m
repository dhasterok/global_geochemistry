function [clon,clat] = smallcircle(lon0,lat0,delta,n)

% start by creating n points on a circle at the pole
clon = [-180:360/n:180]';
clat = 90 - delta*ones(size(clon));

% rotate to a center point at (lon0,lat0)
% it is easier to compute the rotation R from (lon0, lat0) to the pole and
% then use the inverse rotation R' to compute the location of the points
% about (lon0,lat0).
% for the Euler rotation pole, must use a point on the equator and 90deg
% longitude from lon0.  This make rotations simply up or down latitude on a
% constant longitude
R = rotmat(lon0-90,0,90-lat0);

% convert small circle to cartesian
[x,y,z] = sph2cart(clon*pi/180,clat*pi/180,1);
q = [x,y,z]';

% Rotate q
qp = (R'*q)';

% Convert new position of small circle back to spherical coordinates, Q'
[clon,clat] = cart2sph(qp(:,1),qp(:,2),qp(:,3));

clon = clon*180/pi;
clat = clat*180/pi;

return