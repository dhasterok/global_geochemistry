function data = correct_country(data)

addpath worldgrid

% check locations
ind = find(data.latitude < -90 | data.latitude > 90 ...
    | data.longitude < -180 | data.longitude > 180);
if ~isempty(ind)
    warning('LATITUDE or LONGITUDE out of bounds.');
    %ind - this shows what inds its at. Commented out for now.
end

% Dimensions of worldgrid files
nrow = 10801;    % Number of latitudes in worldgrid files
ncol = 21601;   % Number of longitudes in worldgrid files
dl = 180/(nrow - 1);      % Pixel step size in worldgrid files

% Load seafloor age
fprintf('\nLoading oceanic ages...\n');
agefile = 'age.2020.1.GTS2012.1m.nc';
sfage = flipud(ncread(agefile,'z')');

% load elevation
% not currently used
elev = reshape(ncread('ETOPO1_Ice_g_gdal.grd','z'),ncol,nrow)';

fprintf('Computing locations consistent with ages...\n');
% Find index of tectonic age data corresponding to data locations
latind = round((90 - data.latitude)/dl) + 1;
lonind = round((180 + data.longitude)/dl) + 1;

% change indices to points out of bounds to NaN
ind = latind < 1 | latind > nrow | lonind < 1 | lonind > ncol;

data.latitude(ind) = NaN;
data.longitude(ind) = NaN;

latind(ind) = NaN;
lonind(ind) = NaN;

% identify index into sfage
i = sub2ind([nrow ncol],latind,lonind);
j = ~isnan(i);
data.elevation(j) = elev(i(j));

ind = logical(zeros(size(i)));

j = ~isnan(i);
ind(j) = ~isnan(sfage(i(j)));% & elev(i(j)) < 0;
data.country(ind) = {'ocean'};

%figure;
%plot(data.longitude(ind),data.latitude(ind),'.');
if any(strcmp(data.Properties.VariableNames,'avg_age'))
    %ind(j) = ~isnan(sfage(i(j))) & elev(i(j)) < 0 & isnan(data.avg_age(j));
    ind(j) = ~isnan(sfage(i(j))) & isnan(data.avg_age(j));
else
    %ind(j) = ~isnan(sfage(i(j))) & elev(i(j)) < 0;
    ind(j) = ~isnan(sfage(i(j)));
end
%hold on;
%plot(data.longitude(ind),data.latitude(ind),'.');

data.avg_age(ind) = sfage(i(ind));
%figure;
%histogram(sfage(i(ind)))

return
