function [rlon,rlat] = rotsphere(lon,lat,Omega,elon,elat);
% ROTSPHERE - Perfoms spherical rotation.
%
%    [RLON,RLAT] = ROTSPHERE(LON,LAT,OMEGA,ELON,ELAT) rotates the
%    point given by the longitude, LON, and latitude, LAT, by rotation
%    by an angle OMEGA producing the rotated coordinates longitude,
%    RLON, and latitude, RLAT.  LON and LAT must be single numbers or
%    single column or row vectors.
%
% Ref: Cox and Hart, 1986, Plate Tectonics, Blackwell Publishing, p. 227.
% Last Modified: 21 Sept. 2007 by D. Hasterok

% Convert degrees to radians
lon = lon(:)*pi/180;
lat = lat(:)*pi/180;

% Convert input points from spherical to cartesian coordinates
[x,y,z] = sph2cart(lon,lat,1);
A = [x,y,z]';

% Compute spherical rotation matrix
R = rotmat(elon,elat,Omega);

% Compute rotation
Ap = R*A;

% Convert rotated points back to spherical coords
[rlon,rlat] = cart2sph(Ap(1,:),Ap(2,:),Ap(3,:));

% Convert to longitude and latitude in degrees
rlon = rlon(:)*180/pi;
rlat = rlat(:)*180/pi;

return
