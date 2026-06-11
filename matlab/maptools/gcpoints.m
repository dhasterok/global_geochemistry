function [lon,lat,R] = gcpoints(P,Q,n);
% GCPOINTS - Computes a great circle path.
%
%    [LON,LAT] = GCPOINTS(P,Q,N) computes longitude, LON, and
%    latitude, LAT, at N locations along a great circle path
%    between points P and Q.  P and Q must be two element vectors
%    with the first element longitude and second element latitude.
%
%    [LON,LAT,R] = GCPOINTS(P,Q,N) returns the rotation matrix,
%    R, such that P' = R*P where P' is (0,0) and Q' = R*Q' lies
%    on the equator.  Note that the unrotation, P = inv(R)*P' is
%    the same as P = R'*P' since R is orthonormal.
%
% See also: GCPROFILE GCPROJ SPHANGLE
% Last Modified: 26 Sept. 2007 by D. Hasterok

% For testing
%close all;
%clear all;
%
%P = [-114,39];
%Q = [-109,38];
%n = 20;

% All angles are denoted (LON,LAT).

% First: Compute the rotation matricies to move P to (0,0).
%
% This step involves two rotations, the first about the south pole
% to place P on the prime meridan.  The second rotation is an Euler
% pole on the equator.
R1 = rotmat(0,-90,P(1));
R2 = rotmat(90,0,P(2));

% Second: Now compute the new location of Q (Q')
%
% Convert lon/lat to cartesian
[x,y,z] = sph2cart(Q(1)*pi/180,Q(2)*pi/180,1);
q = [x,y,z]';
% Rotate q
qp = R2*R1*q;
% Convert new position of Q back to spherical coordinates, Q'
[lontemp,lattemp] = cart2sph(qp(1),qp(2),qp(3));
Qp = [lontemp lattemp];

% Third: Compute the angle between the equator and Q'.
%
% This equation is the combination of two spherical angle
% computations for a spherical triangle.
%
% Ref: Zwillinger, D. 1996, Standard Mathematical Tables and Formulae,
%      CRC Press, p. 468.
Omega = 180/pi*acos(tan(Qp(1))*cot(acos(cos(Qp(1))*cos(Qp(2)))));
% Seems to be required haven't convinced myself why.
if Qp(2) > 0
    Omega = -Omega;
end

% Fourth: Compute the rotation matrix to bring Q' to the equator
R3 = rotmat(0,0,Omega);

% Single rotation matrix to bring both points onto equator.  Note that
% rotations do not commute (R3*R2*R1 ~= R1*R2*R3).
R = R3*R2*R1;

% Compute the new location of Qp.  Qlat should be 0.
qp = R*q;
[Qlon,Qlat] = cart2sph(qp(1),qp(2),qp(3));
Qp(1) = Qlon;
Qp(2) = Qlat;
if Qlat > 10^-13
    warning('WARNING: Q may not have been moved to the equator properly.');
end

% Compute N points along the equator between P" (0,0) and Q" (lon,0)
Sp = linspace(0,1,n)'*[Qp(1) 0];
% Convert equatorial points to cartesian
[x,y,z] = sph2cart(Sp(:,1),Sp(:,2),1);
sp = [x,y,z]';

% "Unrotate" the points along the equator to the great circle between
% P and Q.
s = (R'*sp)';

% Convert cartesian locations along great circle to spherical
[Slon,Slat] = cart2sph(s(:,1),s(:,2),s(:,3));
lon = Slon*180/pi;
lat = Slat*180/pi;

%Slon = Slon(:);
%Slat = Slat(:);


%figure;
%hold on;
%plot([-180 180],[0 0],'k-');
%plot([0 0],[-90 90],'k-');
%plot(P(1),P(2),'k*',Q(1),Q(2),'k*');
%plot(lon,lat,'r+');
%xlabel('Longitude');
%ylabel('Latitude');
%axis([-180 180 -90 90]);
%set(gca,'Box','on');

return

% Testing stuff
%figure(1);
%hold on;
%plot([-180 180],[0 0],'k-');
%plot([0 0],[-90 90],'k-');
%plot(P(1),P(2),'r*',Q(1),Q(2),'b*');
%[Pxtemp,Pytemp] = rotsphere(P(1),P(2),P(1),0,-90);
%[Qxtemp,Qytemp] = rotsphere(Q(1),Q(2),P(1),0,-90);
%plot(Pxtemp,Pytemp,'r+',Qxtemp,Qytemp,'b+',0,-90,'k^');
%[Pxtemp,Pytemp] = rotsphere(Pxtemp,Pytemp,P(2),90,0);
%[Qxtemp,Qytemp] = rotsphere(Qxtemp,Qytemp,P(2),90,0);
%plot(Pxtemp,Pytemp,'rx',Qxtemp,Qytemp,'bx',90,0,'k^');
%%[x,y,z] = sph2cart(Qxtemp,Qytemp,1);
%%[a,b,c] = sph2cart(Qxtemp,0,1);
%Aytemp = pi/180*Qytemp;
%Axtemp = pi/180*Qxtemp;
%Omega = 180/pi*acos(tan(Axtemp)*cot(acos(cos(Axtemp)*cos(Aytemp))))
%if Aytemp > 0
%    Omega = -Omega;
%end
%[Qxtemp,Qytemp] = rotsphere(Qxtemp,Qytemp,Omega,0,0);
%plot(Qxtemp,Qytemp,'gx');
%xlabel('Longitude');
%ylabel('Latitude');
%axis([-180 180 -90 90]);
%set(gca,'Box','on');
%
%[xp,yp,zp] = sph2cart(pi/180*P(1),pi/180*P(2),1);
%[xq,yq,zq] = sph2cart(pi/180*Q(1),pi/180*Q(2),1);
%sqrt((xp - xq)^2 + (yp - yq)^2 + (zp - zq)^2)
%
%[xp,yp,zp] = sph2cart(pi/180*Pxtemp,pi/180*Pytemp,1);
%[xq,yq,zq] = sph2cart(pi/180*Qxtemp,pi/180*Qytemp,1);
%sqrt((xp - xq)^2 + (yp - yq)^2 + (zp - zq)^2)
%
%break

