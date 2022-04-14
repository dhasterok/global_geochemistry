function P = lognormal(x,mu,sigma)
% LOGNORMAL - Computes a log-normal distribution.
%
%    Computes P(x) for a log-normal distribution P =
%    LOGNORMAL(X,MU,SIGMA) with location parameter, MU, and scale
%    parameter, SIGMA.
%
% Last Modified: 28-July 2011 by D. Hasterok

P = exp(-0.5*(log(x) - mu).^2./sigma.^2)./(x.*sqrt(2*pi)*sigma);

return