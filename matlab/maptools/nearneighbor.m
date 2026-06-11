function [Xi,sd] = nearneighbor(r,R,xi,varargin);
% NEARNEIGHBOR - Nearest Neighbor interpolation.
%
%    XI = NEARNEIGHBOR(r,R,xi) computes the interpolated value,
%    XI given the values, xi, and their distances, r, from a point.
%    The distances to the point are computed prior to calling
%    NEARNEIGHBOR.  The data are weighted by a factor R which is the
%    cutoff distance for interpolation.
%
%    XI = NEARNEIGHBOR(r,R,xi,n) computes interpolated values by
%    computing weights with order n, dependence.
%
%    [XI,SD] = NEARNEIGHBOR(r,R,xi) will also return the standard
%    deviation with the appropriate weights.
%
% Original:  22 Sept. 2007 by D. Hasterok
%
% Corrected standard deviation computation.  The standard deviation was
% previously over estimated.  The standard deviation assumes that all
% measurements have a 20% uncertainty.
% Last Modified: 27 Sept. 2007 by D. Hasterok

if nargin < 3 | nargin > 4
    error('ERROR (nearneighbor.m): Incorrect number of arguments.');
elseif nargin == 4
    n = varargin{1};
else
    n = 2;
end

% Compute weights
w = (1 + 9*r.^n/R^n).^-1;
W = sum(w);

% Compute interpolated point
Xi = sum(xi.*w)/W;

%sd = sqrt(sum((0.2*xi.*w).^2))/sqrt(W);
%sd = sqrt(sum((0.2*xi).^2))/length(xi) + 0.2*mean(xi)*min(r/R);
sd = 0.5*(2 + (min(r)/R)^2)*sqrt(sum((0.2*xi).^2))/length(xi);

return
