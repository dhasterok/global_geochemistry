function rtvert = rtsphere(dlon,dlat)
% vert = rtsphere(dlon,dlat)


dA = rectarea([0 dlon],[0 dlat]);

nlat = 180/dlat;

rtvert = [];
for i = 1:nlat
    latlim = -90 + dlat*[i-1, i];
    ncells = round(rectarea([-180 180],latlim)/dA);

    err(i) = (dA - rectarea([-180 180],latlim)/ncells)/dA;
    dlon_row = 360/ncells;
    temprow = zeros([ncells 4]);
    temprow(:,3) = latlim(1);
    temprow(:,4) = latlim(2);
    if mod(i,2)
        temprow(:,1) = -180 + dlon_row * [0:ncells-1]';
    else
        temprow(:,1) = -180 + dlon_row * [0:ncells-1]' - dlon_row/2;
    end
    temprow(:,2) = temprow(:,1) + dlon_row;
    rtvert = [rtvert; temprow];
end

% errors in the size of the cells relative to 2deg x 2deg at the equator
%err'
%max(abs(err))*100;
%mean(abs(err))*100;

return
