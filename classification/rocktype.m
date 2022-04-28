function ind = rocktype(data,rtype,varargin)

fmt = 'vector';
if nargin == 3
    fmt = varargin{1};
end

switch lower(fmt)
    case 'vector'
        ind = false([height(data) 1]);
    case 'array'
        ind = logical([height(data) length(rtype)]);
    otherwise
        error('Unknown formatting option.');
end

for i = 1:length(rtype)
    switch fmt
        case 'vector'
            ind = ind | strcmp(data.rock_type,rtype{i});
        case 'array'
            ind(:,i) = strcmp(data.rock_type,rtype{i});
    end
end

return