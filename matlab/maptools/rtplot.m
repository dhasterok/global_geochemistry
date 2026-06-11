function rtplot(vert,val)

for i = 1:length(val)
    if isnan(val(i))
        continue;
    end
    f = fill(vert(i,[1 2 2 1 1]),vert(i,[3 3 4 4 3])',val(i));
    f.LineStyle = 'none';
    hold on;
end
p = plotcoast;
p.Color = [0.7 0.7 0.7];

return
