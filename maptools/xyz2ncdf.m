function xyz2ncdf(filename,field,unit,varargin)
% XYZ2NCDF - Converts and XYZ file to NetCDF format.
%
%   xyz2ncdf(filename,field,unit) will convert a file from an X, Y, Value
%   format to NetCDF format.  The outfile will have the same name prefix as
%   the supplied filename with the field names and units supplied as a cell
%   array of strings.  The cell array should not include the longitude or
%   latitude.  The first column should be longitude and the second column
%   latitude.  If latitude and longitude are reversed in the text file,
%   append the function call with string 'Swap', i.e.,
%   xyz2ncdf(filename,field,unit,'Swap').
%
%   Rather than convert a file, the function can write out a data array by
%   adding the option value pairs:
%   xyz2ncdf(filename,field,unit,'Data',x,y,val) where val is a 2-D double
%   array or a cell array of values.e

% Last Modified: 16 Jan. 2022
% D. Hasterok (dhasterok@gmail.com) University of Adelaide

% see if latitude and longitude need to be swapped
swap = false;
opt = 1;
geog = false;
fileflag = true;
while opt + 2 < nargin
    switch lower(varargin{opt})
        case 'swap'
            swap = true;
            opt = opt + 1;
        case 'geographic'
            geog = true;
            opt = opt + 1;
        case 'data'
            fileflag = false;
            x = varargin{opt+1};
            y = varargin{opt+2};
            if iscell(varargin{opt+3})
                v = varargin{opt+3};
            else
                v{1} = varargin{opt+3};
            end
            opt = opt + 4;
        otherwise
            error('Unknown option');
    end
end

if fileflag
    % file for writing
    outfile = [filename(1:end-3),'nc'];
    
    % load data
    data = load(filename);
    
    % assign data to variables, reshape to match the number of latitude and
    % longitude data in the grid.
    n = size(data,2);

    x = unique(data(:,1));
    y = unique(data(:,2));
    for i = 1:n-2
        v{i} = fliplr(reshape(data(:,i+2),length(x),length(y)));
    end
else
    n = 3;
    if ~strcmpi(filename(end-1:end),'nc')
        outfile = [filename,'.nc'];
    else
        outfile = filename;
    end
end

if swap
    yy = x;
    x = y;
    y = yy;

    clear yy;
end

% remove the outfile if it already exists
if exist(outfile,'file')
    fprintf('Removing existing %s.nc first...\n',filename);
    unix(['rm ',outfile]);
end

fprintf('Writing %s.nc...\n',filename);
% create the netCDF file
if geog
    nccreate(outfile, 'lon', 'Dimensions', {'lon',length(x)});
    ncwriteatt(outfile, 'lon', 'standard_name', 'longitude');
    ncwriteatt(outfile, 'lon', 'long_name', 'longitude');
    ncwriteatt(outfile, 'lon', 'units', 'degrees');
    ncwriteatt(outfile, 'lon', '_CoordinateAxisType', 'Lon');
    ncwrite(outfile,'lon',x);

    nccreate(outfile, 'lat', 'Dimensions', {'lat',length(y)});
    ncwriteatt(outfile, 'lat', 'standard_name', 'latitude');
    ncwriteatt(outfile, 'lat', 'long_name', 'latitude');
    ncwriteatt(outfile, 'lat', 'units', 'degrees');
    ncwriteatt(outfile, 'lat', '_CoordinateAxisType', 'Lat');
    ncwrite(outfile,'lat',y);
else
    nccreate(outfile, 'X', 'Dimensions', {'X',length(x)});
    ncwriteatt(outfile, 'X', 'standard_name', 'X');
    ncwriteatt(outfile, 'X', 'long_name', 'X');
    ncwriteatt(outfile, 'X', 'units', 'm');
    ncwriteatt(outfile, 'X', '_CoordinateAxisType', 'X');
    ncwrite(outfile,'X',x);

    nccreate(outfile, 'Y', 'Dimensions', {'Y',length(y)});
    ncwriteatt(outfile, 'Y', 'standard_name', 'Y');
    ncwriteatt(outfile, 'Y', 'long_name', 'Y');
    ncwriteatt(outfile, 'Y', 'units', 'degrees');
    ncwriteatt(outfile, 'Y', '_CoordinateAxisType', 'Y');
    ncwrite(outfile,'Y',y);
end

for i = 1:n-2
    if geog
        nccreate(outfile, field{i+2}, 'Dimensions', {'lon',length(x),'lat',length(y)});
    else
        nccreate(outfile, field{i+2}, 'Dimensions', {'X',length(x),'Y',length(y)});
    end
    ncwriteatt(outfile, field{i+2}, 'standard_name', field{i+2});
    ncwriteatt(outfile, field{i+2}, 'units', unit{i+2});
    ncwrite(outfile, field{i+2}, v{i});
end

return