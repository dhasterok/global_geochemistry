function [dm,F] = erfmodel(m,x)
   
dm = forward(m,x);
F = frechet(m,x);

return

function dm = forward(m,x)
    dm = m(1)*erf(x/m(2));
return

function F = frechet(m,x)
    F = [erf(x/m(2)) ...
         m(1)*2/sqrt(pi)*exp(-(x/m(2)).^2)];
return