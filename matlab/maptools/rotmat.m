function R = rotmat(elon,elat,Omega);
% ROTMAT - Computes a spherical rotation matrix.
%
%    R = rotmat(Elon,Elat,OMEGA) computes the spherical rotation matrix,
%    R, from the Euler pole, E, and the rotation angle, OMEGA.
%
% Last Modified: 21 Sept. 2007 by D. Hasterok

% Convert degrees to radians
elon = elon*pi/180;
elat = elat*pi/180;

Omega = Omega*pi/180;

% Convert Euler pole from spherical to cartesian coordinates
[x,y,z] = sph2cart(elon,elat,1);
E = [x,y,z];

% Rotation Matrix (from Cox and Hart, 1986, Plate Tectonics, p. 227)
R(1,1) = E(1)*E(1)*(1 - cos(Omega)) + cos(Omega);
R(1,2) = E(1)*E(2)*(1 - cos(Omega)) - E(3)*sin(Omega);
R(1,3) = E(1)*E(3)*(1 - cos(Omega)) + E(2)*sin(Omega);
R(2,1) = E(2)*E(1)*(1 - cos(Omega)) + E(3)*sin(Omega);
R(2,2) = E(2)*E(2)*(1 - cos(Omega)) + cos(Omega);
R(2,3) = E(2)*E(3)*(1 - cos(Omega)) - E(1)*sin(Omega);
R(3,1) = E(3)*E(1)*(1 - cos(Omega)) - E(2)*sin(Omega);
R(3,2) = E(3)*E(2)*(1 - cos(Omega)) + E(1)*sin(Omega);
R(3,3) = E(3)*E(3)*(1 - cos(Omega)) + cos(Omega);

return
