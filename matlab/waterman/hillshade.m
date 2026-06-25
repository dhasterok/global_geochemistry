function varargout = hillshade(varargin)
% HILLSHADE - Adds hillshade transparency to an image
%
%   HILLSHADE(x,y,z) creates a hillshade image of a matrix.
%
%   HILLSHADE(...,'Azimuth',angle,'Altitude',angle,'Intensity',factor) will
%   optionally include azimuth, altitude and intensity.
%
%   HILLSHADE(x,y,z,'Histogram',true) will produce a separate figure
%   displaying the distribution of the image gradients.
%
%   hs = HILLSHADE(x,y,z,...) will return the hillshade values.

% parse inputs, they are passed straight through to OctahedralProjection
p = inputParser;
addRequired(p,'X',@isnumeric);
addRequired(p,'Y',@isnumeric);
addRequired(p,'Z',@isnumeric);
addParameter(p,'Azimuth',315,@isnumeric);
addParameter(p,'Altitude',45,@isnumeric);
addParameter(p,'Intensity',1,@isnumeric);
addParameter(p,'Histogram',false,@logical);

parse(p,varargin{:});

histflag = p.Results.Histogram;

% required parameters
x = p.Results.X;
y = p.Results.Y;
z = p.Results.Z;

% lighting azimuth
azimuth = 360 - p.Results.Azimuth + 90;
azimuth(azimuth >= 360) = azimuth - 360;
azimuth = azimuth * (pi/180); %  convert to radians

% lighting altitude
altitude = (90 - p.Results.Altitude)*pi/180; % convert to zenith angle in radians

% lighting intensity
zf = p.Results.Intensity;

clear p

% compute gradient of image
[gmag,gdir] = imgradient(z);

% Histograms 
if histflag
    figure;
    subplot(121);
    histogram(gmag);
    xlabel('Gradient Magnitude');

    subplot(122);
    histogram(gdir);
    xlabel('Gradient Direction');
end

% convert direction to radians
gdir = gdir*pi/180;

% convert magnitude to degrees
grad = atan(zf*gmag);

% ESRIs algorithm
hs = ( (cos(altitude).*cos(grad) ) + ( sin(altitude).*sin(grad).*cos(azimuth-gdir)) );
hs(hs<0) = 0; % set hillshade values to min of 0.

% plot hillshade
figure;
im = imagesc(x,y,z);
im.AlphaData = hs;
prepmap;

if nargout == 1
    varargout{1} = hs;
end

return