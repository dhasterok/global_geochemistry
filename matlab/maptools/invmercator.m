function [x,y] = invmerator(x,y,geoid)


lon = x/Re * 180/pi;
alat = (2*atan(y,Re) - pi/2)*180/pi;
getdatum

return