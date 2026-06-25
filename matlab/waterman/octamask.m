function octamask(varargin)

% parse inputs
p = inputParser;
addParameter(p,'Projection','Waterman',@ischar);
addParameter(p,'Style','Butterfly',@ischar);
addParameter(p,'Inset',false,@islogical);

parse(p,varargin{:});

switch p.Results.Style
    case 'Butterfly'
        if p.Results.Inset
            coords = load('waterman_butterfly_mask_inset.csv');
        else
            coords = load('waterman_butterfly_mask_no_inset.csv');
        end
    case 'M'
    otherwise
        error('Unknown map style.');
end
xm = coords(:,1);
ym = coords(:,2);

fill(xm,ym,'w','EdgeColor','none');
fill(-xm,ym,'w','EdgeColor','none');

return