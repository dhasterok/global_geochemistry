function m = midpt(x,varargin);

% MIDPT - Computes midpoints vector.
%
%    Computes the midpoints of between successive values
%    of a vectors entries.  M = MIDPT(X).
%
% Last Modified: 17 Aug. 2010 by D. Hasterok
%
[nr,nc] = size(x);

dim = 1;
if nargin == 2
    dim = varargin{1};
end

if dim == 2 | nr == 1
    x = x';
end

m = 0.5*(x(1:end-1,:) + x(2:end,:));

if dim == 2 | nr == 1
    m = m';
end

return
