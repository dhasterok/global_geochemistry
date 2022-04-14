function varargout = colormap2(cmap,varargin)
% COLORMAP2 - Additional colormaps
%
%   colormap2(map) will produce a color map either given as a series of
%   color triplets or by given name.  This code adds the following
%   additional colormaps.
%
%   'red', 'blue','green', 'violet', 'brown', 'rwb' (red-white-blue),
%   'ryb' (red-yellow-blue)
if isnumeric(cmap)
    map = palatte(cmap(1,:),cmap(2,:),cmap(3,:),64);
else
    switch cmap
        case 'red'
            map = palatte(0.9*[1 1 1],[1 0.5 0],[0.8627 0.0784 0.2353],64);
        case 'red2'
            map = palatte([1 1 1],[1 0.5 0],[0.8627 0.0784 0.2353],64);
        case 'blue'
            map = palatte(0.9*[1 1 1],[0 0.9804 0.6039],[0 0 0.8039],64);
        case 'green'
            map = palatte(0.9*[1 1 1],[1 0.7843 0.1176],[0.1333 0.5451 0.1333],64);
        case 'violet'
            map = palatte(0.9*[1 1 1],[1 0.6275 0.4784],[0.5804 0 0.8275],64);
        case 'brown'
            map = palatte(0.9*[1 1 1],[0.8706 0.7216 0.5294],[0.5451 0.2706 0.0745],64);
        case 'rwb'
            map = palatte([0.8039 0 0],[1 1 1],[0 0 0.8039],64);
        case 'rwb2'
            map = palatte([0.8039 0.3608 0.3608],[1 1 1],[0.3922 0.5843 0.9294],64);
        case 'ryb'
            map = palatte([0.8039 0 0],[1 0.7843 0.1176],[0 0 0.8039],64);
        case 'ryb2'
            map = palatte([0.8039 0.3608 0.3608],[0.9412 0.9020 0.5490],[0.3922 0.5843 0.9294],64);
        otherwise
            map = cmap;
    end
end

if nargin > 1
    if strcmp(varargin{1},'invert')
        map = flipud(map);
    end
end

if nargout == 0
    colormap(map);
else
    varargout{1} = map;
end

return