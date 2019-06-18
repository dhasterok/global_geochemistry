function [X,Y] = layer(x,y,varargin);

% LAYER   Creates vectors for layered plots.
%    LAYER(X,Y) Creates vectors used to plot constant value
%    layers.  Lengths of X and Y must be equal or differ by one.  Default X 
%    represents top of layer interfaces.  Y represents values
%    within layers.  Last value is value of half-space below
%    bottom interface.
%
%    [X,Y] = LAYER(x,y)  returns vectors X and Y to use with
%    plot.
%
%    [X,Y] = LAYER(x,y,TYPE) allows one to choose whether X
%    represents the top 't' or bottom 'b' of the layers.
%
%    If you wish to plot gradient layers, suggest making a
%    dummy vector for Y.
%
% Modified 4-December-08 by D. Hasterok

% Original 3-April-03 by D. Hasterok
%   4-Dec 2008 Added capability for first element top or bottom
%   of layer

% set inputs
switch nargin
case 2
    type = 'b';
case 3
    type = varargin{1};
otherwise
    error('(layer.m) Incorrect number of arguments.');
end

% make both input vectors columns
x = x(:);
y = y(:);

% length of input vectors
lx = length(x);
ly = length(y);

if abs(lx-ly) > 1
    error('ERROR: Vectors X and Y must be same length or different by 1.');
    return
end

if lx == ly + 1
    X = [x(1:lx-1)'; x(2:lx)'];
    X = X(:);

    Y = [y'; y'];
    Y = Y(:);

    return;

elseif lx + 1 == ly
    X = [x'; x'];
    X = X(:);

    Y = [y(1:ly-1)'; y(2:ly)'];
    Y = Y(:);

    return;
end

X = zeros(2*lx,1);
Y = zeros(2*ly,1);

switch lower(type);
    case 'b'
        X(1:2:2*lx-1) = x(1:lx);
        X(2:2:2*lx-2) = x(2:lx);
        % add bottom of bottom layer
        X(2*lx) = x(lx) + 0.1*(x(lx) - x(1));
    case 't'
        % add top of top layer
        X(1) = x(1) - 0.5*(x(2) - x(1));
        % check to make sure top layer does not cross zero
        if X(1) < 0
            X(1) = 0;
        end
        X(2:2:2*lx) = x(1:lx);
        X(3:2:2*lx-1) = x(1:lx-1);
    otherwise
        error('Incorrect type of layering');
end

Y(1:2:2*ly-1) = y(1:ly);
Y(2:2:2*ly)   = y(1:ly);

return
