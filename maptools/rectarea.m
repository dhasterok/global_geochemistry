function A = rectarea(lonlim,latlim)
% A = rectarea(lonlim,latlim)


R = 6371;

lat1 = deg2rad(min(latlim));
lat2 = deg2rad(max(latlim));
lon1 = deg2rad(min(lonlim));
lon2 = deg2rad(max(lonlim));

A = 2*pi*R^2 * abs(sin(lat1) - sin(lat2)) .* abs(lon1 - lon2)/360;

return
