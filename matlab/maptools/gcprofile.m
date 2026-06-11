function [plon,plat,Xi,sd,pts] = gcprofile(lon,lat,xi,R,N,P1,P2,varargin);
% GCPROFILE - Computes values along a great circle.
%
%    [plon,plat,Xi] = GCPROFILE(LON,LAT,XI,R,N,P1,P2) computes the
%    interpolated data given by the vector defined by the triplet
%    longitude, latitude and magnitude (LON,LAT,XI).  The interpolated
%    values, Xi, are computed using a nearest neighbor scheme, using all
%    points that fall within a radius, R, of the great circle defined by
%    the points P1 and P2.  The program interpolates for N points along
%    the great circle located at (PLON, PLAT).
%
%    Note: The radius, R, should be normalized to the radius of the
%    sphere.  Thus, for the Earth, the distance from the line to the
%    points should be given as R = r/6371.
%
%    [plon,plat,Xi] = GCPROFILE(LON,LAT,XI,R,N,P1,P2,ORDER) will compute
%    interpolated values by calling NEARNEIGHBOR with a weighting
%    function to order, ORDER, dependence.
%
% See also: GCPOINTS NEARNEIGHBOR
% Last Modified: 22 Sept. 2007 by D. Hasterok

% Compute N points on the great circle path from points P1 to P2
[plon,plat] = gcpoints(P1,P2,N);

if nargin < 7 | nargin > 8
    error('ERROR (gcprofile.m): Invalid number of arguments.');
elseif nargin == 4
    order = varargin{1};
else
    order = 2;
end

d2r = pi/180;

% Convert points to cartesian coordinates
[x,y,z] = sph2cart(d2r*lon(:),d2r*lat(:),1);

% Convert point along line to cartesian
[xp,yp,zp] = sph2cart(d2r*plon,d2r*plat,1);

Xi = zeros(size(plon));
sd = zeros(size(plon));
for i = 1:N
    % Compute distance between the data and line
    r = sqrt((x - xp(i)).^2 + (y - yp(i)).^2 + (z - zp(i)).^2);

    % Pick all points <= R from the data 
    ind = find(r <= R);
    pts{i} = ind;

    % Perfom nearest-neighbor interpolation scheme
    if isempty(ind)
        Xi(i) = NaN;
        sd(i) = NaN;
        continue
    end
    [Xi(i),sd(i)] = nearneighbor(r(ind),R,xi(ind),order);
end

return
