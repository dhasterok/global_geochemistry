function [clon,clat,ang] = gcproj(lon,lat,P,Q);
% GCPROJ - Projects a set of points onto a great circle.
%
%    [GCLON,GCLAT] = GCPROJ(LON,LAT,P,Q) projects a set of points
%    (LON,LAT) onto a great circle defined by P and Q.  GCPROJ
%    returns the longitude, GCLON, and latitude, GCLAT, of the
%    original points at the intersection of a great circle
%    orthogonal to the reference great circle.
%
%    [GCLON,GCLAT,ANGDIST] = GCPROJ(LON,LAT,P,Q) will return
%    the angular distance, ANGDIST (in degrees), of each point to
%    the great circle in addition to the intersection point.  
%
% See also: GCPOINTS GCPROFILE
% Last Modified: 26 Sept. 2007 by D. Hasterok

% Compute rotation matrix to bring great circle to equator
% and all other points relative to the equator
[dummy(:,1),dummy(:,2),R] = gcpoints(P,Q,2);

% Convert from degrees to radians
lamda = pi/180*lon(:);
phi = pi/180*lat(:);

% Convert points from spherical to cartesian
[X(:,1),X(:,2),X(:,3)] = sph2cart(lamda,phi,1);

% Rotate points, setting them relative to the equator
Xp = R*X';

% Since the x- and y-coordinates define the latitude and the z-coordinate
% defines the longitude.  Setting the z-coordinate to zero will shift the
% points to the equator.  This can be done to "project" the points onto the
% great circle because the great circle defined by the longiude is always
% orthogonal to the equator.

% The long way:
%[plambda,pphi] = cart2sph(Xp(1,:),Xp(2,:),Xp(3,:));
%plon = 180/pi*plambda';
%plat = 180/pi*pphi';
%
%[x,y,z] = sph2cart(plambda,zeros(size(plambda)),1);
%Xeq = [x,y,z];

% The short way:
Xeq = [Xp(1:2,:); 0*Xp(3,:)];

% The latitude in the equatorial reference frame represents the angular
% distance to the great circle.  The latitude can be simply computed from
% the arccosine of the z-coordinate since the distance from the center is
% unity.
ang = 180/pi*asin(Xp(3,:))';

% Unproject points to the original great circle
Xun = R'*Xeq;

% Convert from cartesian to spherical
[clambda,cphi] = cart2sph(Xun(1,:),Xun(2,:),Xun(3,:));

clon = 180/pi*clambda';
clat = 180/pi*cphi';

return
