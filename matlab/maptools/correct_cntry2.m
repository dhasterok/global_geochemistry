function data = correct_country(data)

% check locations
ind = find(data.latitude < -90 | data.latitude > 90 ...
    | data.longitude < -180 | data.longitude > 180);
if ~isempty(ind)
    warning('LATITUDE or LONGITUDE out of bounds.');
    %ind - this shows what inds its at. Commented out for now.
end

% Dimensions of worldgrid files
dl = 1/30;      % Pixel step size in worldgrid files
nrow = 5401;    % Number of latitudes in worldgrid files
ncol = 10801;   % Number of longitudes in worldgrid files

% Load tectonic age
fprintf('\nLoading oceanic ages...\n');
sfage = worldgrid('sfage');
[iage,ind] = addages;
sfage(ind) = iage;

fprintf('\nComputing locations consistent with ages...\n');
% Find index of tectonic age data corresponding to data locations
latind = round((90 - data.latitude)/dl) + 1;
lonind = round((180 + data.longitude)/dl) + 1;

ind = ~isnan(sfage(latind,lonind));
data.country{ind} = 'ocean';

ind = isnan(data.avg_age);
data.avg_age{ind} = sfage(latind(ind),lonind(ind));

return
