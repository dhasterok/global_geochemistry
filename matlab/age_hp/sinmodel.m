function [dm,F] = sinmodel(m,x)
   
dm = forward(m,x);
F = frechet(m,x);

return
 
function dm = forward(m,x)
    dm = m(1)*sin(2*pi*x/m(2) + m(3)) - m(4);
return

function F = frechet(m,x)
    F = [sin(2*pi*x/m(2) + m(3)) ...
        -2*pi*m(1)*m(2)^-2*x.*cos(2*pi*x/m(2) + m(3)) ...
        2*pi*m(1)*cos(2*pi*x/m(2) + m(3)) ...
        -ones(size(x))];
return