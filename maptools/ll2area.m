function area = ll2area(glat,dl);

% determin authalic latitude and Earth's radius
[alat,Re] = authalic([glat-dl/2 glat+dl/2],'nad83');
Re = Re/1e3;

% convert to radians
alat = alat*pi/180;

n = 360/dl;

area = 2*pi*Re^2*(sin(alat(:,2)) - sin(alat(:,1)))/n;

return