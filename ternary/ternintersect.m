function [ai,bi,ci] = ternintersect(a1,b1,c1, a2,b2,c2)
% TERNINTERSECT - determines interesction of two lines on a ternary plot
%
%   [ai,bi,ci] = ternintersect(a1,b1,c1, a2,b2,c2) computes the
%   intersection coordinate triplet (ai,bi,ci) given two lines defined by
%   the ternary endpoints (a1,b1,c1) and (a2,b2,c2).
%
%   Note the returned values are normalized [0,1].
%
% D. Hasterok, 18 Nov 2021

s = a1(1) + b1(1) + c1(1);

% convert the endpoints of ternary lines to cartesian points
[x1,y1] = tern2xy(a1,b1,c1);
[x2,y2] = tern2xy(a2,b2,c2);

% line 1 parameters
m11 = diff(y1)/diff(x1);    % slope
m12 = -m11*x1(1) + y1(1);   % intercept

% line 2 parameters
m21 = diff(y2)/diff(x2);    % slope
m22 = -m21*x2(1) + y2(1);   % intercept

% intersection point
xi = (m22 - m12)/(m11 - m21);
yi = m11*xi + m12;

% convert intersection point from cartesian to ternary
[ai,bi,ci] = xy2tern(xi,yi);

ai = ai*s;
bi = bi*s;
ci = ci*s;

return