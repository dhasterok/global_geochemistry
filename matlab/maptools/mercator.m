function [x,y] = mercator(lat,lon,geoid)
% MERCATOR - computes mercator coordinates
%
%   [x,y] = mercator(lat,lon,geoid) performs a Mercator projection for the
%   latitude and longitude locations, correcting for the given geoid.
%
% see also: invmercator

% compute authalic latitude and Earth radius
[alat,Re] = authalic(lat,geoid);

% compute mercator coordinates
x = Re*lon*pi/180;
y = Re*log(tan(pi/4 + alat*pi/360));

return