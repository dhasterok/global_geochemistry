function k2o = getK2O(data)

if ~any(strcmp('k_ppm',fieldnames(data)))
    k2o = data.k2o;
    return;
end

k2o = nan([length(data.sio2),1]);
for i = 1:length(data.sio2)
    if isnan(data.k_ppm(i))
        k2o(i) = data.k2o(i);
    elseif isnan(data.k2o(i))
        k2o(i) = (2*39.0986 + 15.9994)/(2*39.0986*1e4)*data.k_ppm(i);
    elseif ~isnan(data.k2o(i)) && ~isnan(data.k_ppm(i))
        if data.k2o(i) < 1
            k2o(i) = (2*39.0986 + 15.9994)/(2*39.0986*1e4)*data.k_ppm(i);
        else
            k2o(i) = data.k2o(i);
        end
    end
end

return
