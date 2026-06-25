function fillhex(xv,yv,h)

zbottom = zeros(1,7);
ztop = h * ones(1,7);

fill3(xv,yv,zbottom,h);
fill3(xv,yv,ztop,h);

for i = 1:6
    fill3([xv(i), xv(i+1), xv(i+1), xv(i)], ...
        [yv(i), yv(i+1), yv(i+1), yv(i)], ...
        [zbottom(i), zbottom(i), ztop(i+1), ztop(i+1)],h);
end

return