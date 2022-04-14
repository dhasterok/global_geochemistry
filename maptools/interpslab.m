function vdepth = interpslab(slab,vlon,vlat)

lon = [];
lat = [];
depth = [];
for i = 1:length(slab)
    for j = 1:length(slab(i).contour)
        pos = slab(i).contour(j);
        ind = pos.lon > 180;
        pos.lon(ind) = pos.lon(ind)-360;
        
        lon = [lon; pos.lon];
        lat = [lat; pos.lat];
        depth = [depth; repmat(slab(i).depth(j),size(pos.lon))];
    end
end

F = scatteredInterpolant(lon,lat,depth,'linear','none');

vdepth = F(vlon,vlat);

return