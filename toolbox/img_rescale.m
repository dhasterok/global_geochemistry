function [X,Y] = img_rescale(fig,nx,ny);
% IMG_RESCALE - Rescales an image.
%
%   IMG_RESCALE maps the coordinates of an image to a new set of axes.
%
%   Example call:
%       [X,Y] = IMG_RESCALE(FIG,NX,NY)
%
%   Inputs:
%       FIG - Figure number or name
%       NX,NY - Resolution of the image in the X and Y dimensions
%
%   Outputs:
%       X,Y - Coordinates of the rescaled image
%
% Last Modified: 21-Feb 2007
% See also: FIGDIG

% Scale the image
fprintf('Select control points.\n\n');
[x,y] = getpoints(fig);
if length(x) < 2
    error('ERROR: Not enough points to determine scale.');
end

if length(x) == 2
    if x(1) == x(2) | y(1) == y(2)
        error('ERROR: Must choose points with different x''s or y''s.');
    end
end

str = sprintf('Control points.\n  Enter each point pair in order separated by a '';''.  Example: [0 0; 1 1; -1 2]\n\n  cp = ');
%cp = input('Enter control points: ');
cp = input(str);
xp = cp(:,1);
yp = cp(:,2);

if length(xp) ~= length(x)
    error('ERROR: Incorrect number of points entered.');
end

% Find the linear scaling factors
A = [ones(size(x)) x];
mx = inv(A'*A)*A'*xp;
B = [ones(size(y)) y];
my = inv(B'*B)*B'*yp;

% Rescale the image
X = mx(1) + mx(2)*[0:nx];
Y = my(1) + my(2)*[0:ny];

return
