function varargout = lngammamodel(m,x)
% fits a log gamma PDF given the positions, x, and scale parameters, 
% alpha = m(1) and beta = m(2).

varargout{1} = forward(m,x);
if nargout == 2
    varargout{2} = frechet(m,x);
end

return

function dm = forward(m,x)

dm = lngamma(x,m(1),m(2));

return

function F = frechet(m,x)

arg = m(2)*x - exp(x)/m(1);
F = [m(1)^(-(m(2) + 2))*exp(x + arg) - m(2)*m(1)^(-(m(2)-1)), ...
    m(1)^-m(2)*exp(arg).*(-log(m(1)) - psi(0,m(2)) + x)]/gamma(m(2));

return