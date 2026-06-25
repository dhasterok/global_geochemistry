function c = colorbyvalue(val,bounds,cmap)

% nomalize dataset to bounds
val = (val - bounds(1))/(bounds(2) - bounds(1));
val(val < 0) = 0;
val(val > 1) = 1;

c = zeros(height(val),3);
for i = 2:size(cmap,1)
    ind = cmap(i-1,1) <= val & val <= cmap(i,1);

    x = (val(ind) - cmap(i-1,1));
    for j = 1:3
        c(ind,j) = (cmap(i,j+1) - cmap(i-1,j+1))/(cmap(i,1) - cmap(i-1,1))*x + cmap(i-1,j+1);
    end
end

return