function cntry = correct_cntry(filename);

fid = fopen(filename);
places = textscan(fid,'%f%f%f%s\n','Delimiter',',');

id = places{1};
lat = places{2};
lon = places{3};
cntry = places{4};

% check locations
ind = find(lat < -90 | lat > 90 ...
    | lon < -180 | lon > 180);
if ~isempty(ind)
    warning('LATITUDE or LONGITUDE out of bounds.');
    ind
end

% Dimensions of worldgrid files
dl = 1/30;      % Pixel step size in worldgrid files
nrow = 5401;    % Number of latitudes in worldgrid files
ncol = 10801;   % Number of longitudes in worldgrid files

% Load tectonic age
%fprintf('\nLoading continent mask...\n');
%oc = worldgrid('contmask');
%[tage,alat,alon] = worldgrid('tage');
%[alon,alat] = meshgrid(alon,alat);

age = worldgrid('sfage');
[iage,ind] = addages;
age(ind) = iage;

% Find index of tectonic age data corresponding to heat flow locations
latind = round((90 - lat)/dl) + 1;
lonind = round((180 + lon)/dl) + 1;

for i = 1:length(id)
%    if oc(latind(i),lonind(i))
%        cntry{i} = 'ocean';
%    end
    if ~isnan(age(latind(i),lonind(i)))
        cntry{i} = 'ocean';
    end
end

ind = strcmp('ocean',cntry);
figure;
plotcoast;
plot(lon,lat,'b.');
plot(lon(ind),lat(ind),'r.');

fid = fopen([filename(1:end-4),'oc.csv'],'w');
for i = 1:length(id)
    fprintf(fid,'%i,%f,%f,%s\n',id(i),lat(i),lon(i),cntry{i});
end
fclose(fid); 

return
