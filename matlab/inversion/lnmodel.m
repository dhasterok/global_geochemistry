% --------------------------------------------------------------------
% you may change the function name, but not the inputs/outputs
function [varargout] = model(m,x);

% --------------------------------------------------------------------
% Do not edit the code between this and the next line below

% compute model
dm = forward(m,x);
varargout{1} = dm;

% only compute frechet and weighting matricies if needed
if nargout > 1 | nargout < 4
    % compute frechet matrix
    F = frechet(m,x,dm);
    varargout{2} = F;

    if nargout == 3
        % compute data and model weighting
        Km = weighting(m,x,dm);
        varargout{3} = Km;
    end
elseif nargout ~= 1
    error('(MODEL): Incorrect number of output arguments requested.');
end

return

% --------------------------------------------------------------------

% You may change the functions below this point, but not the function
% names as they are used above.


% FORWARD - computes synthetic data/forward model based on a prescribed
% (empirical model or physics)
%
%  d_m = f(x1,x2,...) = forward(m,x)
%
%  where each column in x is a different independent variable and m is
%  the different model constants/coefficients to be determined by
%  inversion.
function dm = forward(m,x);
    dm = exp( - (log(x) - m(1)).^2 ./ (2*m(2)^2) ) ./ ( x * sqrt(2*pi) * m(2) );
return


% FRECHET - aka Jacobian (matrix of derivatives w.r.t. model parameters)
%
%      [ df    df      ]
%  F = [ ---   --- ... ]
%      [ dm1   dm2     ]
%
%  where each column represents the derivative of the forward model
%  w.r.t. a particular model (fitting) parameter.  Each row is
%  evaluated w.r.t. a different combination of input parameters.
function F = frechet(m,x,dm);
    a = (log(x) - m(1)) / m(2);

    F = y.*[ a / m(2) , (a.^2 - 1) / m(2) ];
return


% WEIGHTING - used to determine weight misfit fuctional by uncertainties
% on the independent variables
%
%       [ df    df      ]
%  Kx = [ ---   --- ... ]
%       [ dx1   dx2     ]
%
%  where each column represents the derivative of the forward model
%  w.r.t. a particular independent (known input) variable.  Each row is
%  evaluated w.r.t. a different combination of input parameters.
function Kx = weighting(m,x,dm);
    Kx = y*(log(x) - m(1) + m(2)^2) ./ (x * m(2)^2);
return

