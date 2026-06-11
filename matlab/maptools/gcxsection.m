function [dx,ind] = gcxsection(lon,lat,r,P,Q)
% GCXSECTION - Computes the positions of points along a cross-section
%
%   [dx,ind] = gcxsection(lon,lat,r,P,Q) computes the distance of
%   observation points within a distance r of, and projected to a great 
%   circle path from P to Q.
%
% 24 Feb 2022 by D. Hasterok

C = Constants;

[plon,plat,angle] = gcproj(lon,lat,P,Q);

dx = C.Rearth*sphangle(plon,plat,P(1),P(2));
dq = C.Rearth*sphangle(Q(1),Q(2),P(1),P(2));

ind = abs(angle*pi/180*C.Rearth) <= r & 0 <= dx & dx <= dq;

return