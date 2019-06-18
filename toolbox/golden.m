function golden(varargin);
% GOLDEN - plots golden ratio aspect ratio.
%
% GOLDEN plots within a box with aspect ratio with x:y = 1.618:1.
%
% GOLDEN(AX) plots within a box with the axis (X, Y or Z) as 1.618 and all
% others relative length 1.

ar = [1 1 1];
ax = 'x';

% golden ratio
phi = 1.618;

if nargin == 1
    ax = varargin{1};
    if ~strcmp(ax,'x') & ~strcmp(ax,'y') & ~strcmp(ax,'z')
        warning('Unknown axis.');
        ax = 'x';
    end
end

switch ax
    case 'x'
        ar(1) = phi;
    case 'y'
        ar(2) = phi;
    case 'z'
        ar(3) = phi;
end

pbaspect(ar);

return
