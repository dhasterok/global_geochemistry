function P = lognormal(x,alpha,beta)
% LNGAMMA - Computes a log-gamma distribution.
%
%    Computes P(x) for a log-gamma distribution P =
%    LNGAMMA(X,alpha,beta) with location parameter, alpha, and scale
%    parameter, beta.
%
% Last Modified: 22-September 2020 by D. Hasterok

P = (exp(beta*x).*exp(-exp(x)/alpha))./(alpha^beta*gamma(beta));

return