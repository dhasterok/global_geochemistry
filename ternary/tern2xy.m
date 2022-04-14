function [x,y] = tern2xy(a,b,c);
% TERN2XY - Converts ternary points to cartesian
%
%   [x,y] = tern2xy(a,b,c) converts ternary points (a,b,c) to cartesian
%   (x,y) coordinates.  The cartesian origin is defined as the midpoint on
%   the B-C edge.
%
%   The order of the axes are as follows:
%
%                    A
%                   / \
%                  /   \
%                 B -•- C
%
%   where the cartesian origin defined at '•' above.
%
%   See also XY2TERN

% half-width
w = 0.5;

% vertical scale
h = 0.5/tan(pi/6);

% nomalize (a,b,c)
s = (a + b + c);
a = a./s;
b = b./s;
c = c./s;

% convert to xy
y = a*h;
x = (1 - b)*h/cos(pi/6) - y*tan(pi/6) - w ;

return
