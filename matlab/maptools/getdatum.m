function [a,b] = getdatum(datum)
% GETDATUM - Gets constants for map ellipsoids.
%
%    [A,B] = GETDATUM(DATUM) gets the semimajor, A, and semiminor, B,
%    axes for DATUM (string).
%
%    Datums:
%       nad27, nad83, wgs84, grs80, wgrs80, wgrs84, clarke1866,
%       clarke1880, australia1965, australia1984, krasovsky1940,
%       intl1924, hayford1909, airy1830, bessel1841, everest1830,
%       nzmg, nztm
%
% Last Modified: 19 Sept. 2007 by D. Hasterok
datum = lower(datum);

switch datum
    case {'nad27','clarke1866'}
        a = 6378206.4;
        b = 6356583.8;
    case {'nad83','wgs84','grs80','wgrs80','wgrs84','nztm'}
        a = 6378137.0;
        b = 6356752.3;
    case 'australian1965'
        a = 6378160.0;
        b = 6356774.7;
    case {'australian1984'; 'gda84'}
        a = 6378160.0;
        b = 6356774.7;
    case {'australian1994'; 'gda94'}
        a = 6377137.0;
        b = 6355755.7;
    case 'krasovsky1940'
        a = 6378245.0;
        b = 6356863.0;
    case {'intl1924','hayford1909','nzmg'}
        a = 6378388.0;
        b = 6356911.9;
    case 'clarke1880'
        a = 6378249.1;
        b = 6356514.9;
    case 'airy1830'
        a = 6377563.4;
        b = 6356356.9;
    case 'bessel1841'
        a = 6377397.155;
        b = 6356079.0;
    case 'everest1830'
        a = 6377376.3;
        b = 6356075.4;
    otherwise
        error('ERROR: Datum not recognized.');
end

return
