function [ecc,a,b] = eccentricity(varargin);
% ECCENTRICITY - Computes the eccentricity of an ellipse.
%
%    E = ECCENTRICTY(A,B) computes the eccentricity E of an elipse
%    defined by a semimajor axis A and semiminor axis B.
%
%    E = ECCENTRICTY(DATUM) computes the eccentricity E of the ellipsoid
%    DATUM.
 
% Change log:
%    Modified: 23 Oct. 2012 by D. Hasterok, fixed a bug with input
%    values.
%
%    Modified: 3 Apr. 2012 by D. Hasterok, added shortcut to compute
%    eccentricity using dataum
%
%    Modified: 19 Sept. 2007 by D. Hasterok
if nargin == 1
    datum = varargin{1};
    [a,b] = getdatum(datum);
elseif nargin == 2
    a = varargin{1};
    b = varargin{2};
else
    error('ERROR (eccentricity): Incorrect number of arguments.');
end

ecc = sqrt(1 - b^2/a^2);

return
