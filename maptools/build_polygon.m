function arcpoly = build_polygon(idl_flag,vlon,vlat,delta)
% arcpoly = build_polygon(idl_flag,vlon,vlat,delta)
% constructs a polygon from a set of points points

for j = 1:length(vlon)
    [clon,clat] = smallcircle(vlon(j),vlat(j),delta,20);
    if idl_flag
        clon(clon > 90) = clon(clon > 90) - 360;
    end
    ptemp(j) = polyshape(clon,clat);

    if j == 2
        arcpoly = union(ptemp(1),ptemp(2));
    elseif j > 2
        arcpoly = union(arcpoly,ptemp(j));
    end
end

if idl_flag
    tmppoly = arcpoly;
    tmppoly.Vertices(:,1) = 360 + tmppoly.Vertices(:,1);

    arcpoly = union(arcpoly,tmppoly);
end

return