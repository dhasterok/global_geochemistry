

function plotslab(slab)

C = parula(12);
depth = [];
p = [];
slab(end)
c = 1;
for i = 1:length(slab)
    for j = 1:length(slab(i).contour)
        ind = find(depth == slab(i).depth(j));
        if isempty(ind)
            depth = [depth; slab(i).depth(j)];
            ltxt{c} = num2str(depth(c));
            ind = length(depth);
            c = c + 1;
        end
        tmp = plot(slab(i).contour(j).lon,slab(i).contour(j).lat, ...
                '-','Color',C(ind,:));
        hold on;
        if length(p) < length(depth)
            p = [p; tmp];
        end
    end
end

legend(p,ltxt);

return