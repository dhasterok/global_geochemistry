function rt = rtstats(dlon,dlat,lon,lat,data);

rt.vert = rtsphere(dlon,dlat);
n = length(rt.vert(:,1));

rt.n = zeros([n,1]);
rt.median = zeros([n,1]);
rt.mean = zeros([n,1]);
rt.std = zeros([n,1]);
for i = 1:n
    if rt.vert(i,1) < -180
        in = (-180 <= lon & lon < rt.vert(i,2) ...
            & rt.vert(i,3) <= lat & lat < rt.vert(i,4)) | ...
            (360-rt.vert(i,1) <= lon & lon < 180 ...
            & rt.vert(i,3) <= lat & lat < rt.vert(i,4));

    else
        in = rt.vert(i,1) <= lon & lon < rt.vert(i,2) ...
            & rt.vert(i,3) <= lat & lat < rt.vert(i,4);
    end
    rt.n(i) = sum(in);
    if rt.n(i) ~= 0
        rt.median(i) = median(data(in));
        rt.mean(i) = mean(data(in));
        rt.std(i) = std(data(in));
    else
        rt.median(i) = NaN;
        rt.mean(i) = NaN;
        rt.std(i) = NaN;
    end
end

return
