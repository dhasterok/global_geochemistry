function [x,y] = figdig(varargin);
% FIGDIG - Digitizes a figure.
%
%    FIGDIG digitizes a figure or images.  If an image is selected you will
%    be prompted to select a few control points before digitizing.
%
%    Example calls:
%        To digitize an existing figure, use [X,Y] = FIGDIG(FIG). or to 
%        digitize an image file, [X,Y] = FIGDIG(FILE).
%
%    Inputs:
%        FIG  - Figure number/handle
%        FILE - Name of the image file (must be a type readable by
%               MATLAB.  See IMREAD for known types).
%
%    Ouput:
%        X,Y - Digitized coordinates
%
% Last Modified: 21-Feb 2007
% See also IMREAD, IMG_RESCALE, GETPOINTS

if isa(varargin{1},'double')
    % Select figure
    fig = varargin{1};

    figure(fig);
elseif ischar(varargin{1})
    filename = varargin{1};
    len = length(filename);

    % Read the image.  The first argument is the name of the file, the
    % second is the type.  See IMREAD for types.
    try
        im = imread(filename);
    catch
        im = imread(filename,filename(len-2:len));
    end
    fig = figure;

    % Plot the image
    imagesc(im);
    [ny,nx,nz] = size(im);

    % Rescale the image
    [X,Y] = img_rescale(fig,nx,ny);

    % Replot the rescaled image
    fig = figure;
    imagesc(X,Y,im);
    axis xy;
end

fprintf('Digitizing.\n');
[x,y] = getpoints(fig);

return
