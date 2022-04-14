function p = staircase(edges,input_type,data,varargin);
% STAIRCASE - Makes a staircase plot.
%
%   p = staircase(edges,height)
scale = 'linear';
colour = [];
opt = 1;
while opt + 3 < nargin
    switch lower(varargin{opt})
        case 'scale'
            scale = varargin{opt+1};
        case 'color'
            colour = varargin{opt+1};
            if size(colour,1) == 1
                colour = repmat(colour,size(data,2),1);
            end
        case 'displaystyle'
            displaystyle = varargin{opt+1};
        otherwise
            error('Unknown option.');
    end
    opt = opt + 2;
end

switch lower(input_type)
    case 'height'
        height = data;
    case 'values'
        height = zeros(1,length(edges)-1);
        for i = 1:length(edges)-1
            height(i) = sum(edges(i) <= data & data < edges(i+1));
        end
    otherwise
        error('Unknown input_type.');
end

% change height of 
switch lower(scale)
    case 'linear'
        % do nothing
    case 'log'
        height = log10(height);
        ind = isnan(height) | height == -Inf;
        height(ind) = nan;
        height(ind) = -1;
    otherwise
        error('Unknown scale option.');
end     

[nr,nc] = size(height);

if nr == 1 | nc == 1
    height = height(:);
    nc = 1;
elseif nc == length(edges)
    height = height';
    nc = nr;
end
edges = edges(:);

x = [edges edges]';
x = x(:);

for i = 1:nc
    y = [[0; height(:,i)] [height(:,i); 0]]';
    y = y(:);

    if strcmpi('displaystyle','stacked')
        if i == 1
            yy = y;
        else
            yy = yy + y;
        end
        p(i) = plot(x,yy,'-');
    else
        p(i) = plot(x,y,'-');
    end
    if ~isempty(colour)
        set(p(i),'Color',colour(i,:));
    end
end

% create semilogy ticks and tick labels
if strcmpi(scale,'log')
    hpax([floor(min(height,[],'all')) ceil(max(height,[],'all'))],'y');
    ylabel('');
end

return