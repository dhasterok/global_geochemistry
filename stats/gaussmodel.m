function varargout = gaussmodel(m,x)
% fits a gaussian CDF given the positions, x, and mean, m(1), standard
% deviation, m(2).

varargout{1} = forward(m,x);
if nargout == 2
    varargout{2} = frechet(m,x);
end

return

function dm = forward(m,x)

dm = 0.5*(1 + erf((x - m(1))/(sqrt(2)*m(2))));

return

function F = frechet(m,x)

arg = (x - m(1))/(sqrt(2)*m(2));
F = -[1 m(1)/m(2)].*repmat(exp(-arg.^2)/(sqrt(pi)*sqrt(2)*m(2)),1,2);

return