function [map,cb] = gtspal(varargin)

if nargin > 0
    filename = varargin{1};
else
    filename = 'gts_limited.xlsx';
end

ts = readtable(filename);

ts = sortrows(ts,'age_lower');

t = [min(ts.age_lower):1:max(ts.age_upper)];


c = 1;
for i = 1:length(t);
    if ts.age_upper(c) < t(i)
        c = c + 1;
        continue;
    end

    map(i,:) = ts{c,{'R','G','B'}};
end

map = map./255;

colormap(map);
caxis([ts.age_lower(1) ts.age_upper(end)]);
cb = colorbar;
cb.TickDirection = 'out';
cb.Ticks = [ts.age_lower; ts.age_upper(end)];
cb.YLabel.String = 'Age [Ma]';

return
