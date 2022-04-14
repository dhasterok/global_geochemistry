function r = lnrand(mu,sigma);
% LNRAND - Log normally distributed random numbers
%
%   R = LNRAND(MU,SIGMA) produces a random number R with the standard scale
%   factors MU and SIGMA.

r = randn(1);
r = exp(mu + r*sigma);

return
