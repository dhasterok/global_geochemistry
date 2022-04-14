function varargout = gaussian(x,varargin)

% GAUSSIAN - Gaussian (normal) distribution.
%
%    PDF = gaussian(x,mu,sigma) computes the probablility distribution function
%    for a gaussian distribution where mu is the mean and sigma is the standard
%    deviation of the random variable x.
%
%    [PDF,W] = gaussian(x,mu,sigma) additionally computes the Shipiro-Wilk
%    statistic as a test of normality.
%
%    [mu,sigma] = gaussian(x) computes the mean, mu, and sample standard
%    deviation, sigma, given the random variable x.

% Modified 21-Oct 2014 by D. Hasterok to include the Shapiro-Wilk statistic and
% compute mean and standard deviation if they are presently unknown.


if nargin == 3
    mu = varargin{1};
    sigma = varargin{2};

    PDF = exp(-0.5*(x - mu).^2./sigma.^2)./(sqrt(2*pi)*sigma);
    varargout{1} = PDF; 

    if nargout == 2
        W = shapirowilk(x);
        varargout{2} = W;
    elseif nargout > 2
        error('Too many output arguements requested.');
    end
elseif nargin == 1
    n = length(x);

    mu = sum(x) / n;
    sigma = sqrt(sum((x - mu).^2) / (n - 1));

    varargout{1} = mu;
    varargout{2} = sigma;
    
else
    error('Incorrect number of input arguements.');
end

return

