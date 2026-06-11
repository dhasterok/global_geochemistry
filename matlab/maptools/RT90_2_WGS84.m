%RT90_2_WGS84 Convert coordinate from RT90 to WGS84 format.
%
%   WGS84 = RT90_2_WGS84(RT90) 
%
%       RT90          N by 2: coordinates in RT90 (x and y)
%       WGS84         N by 2: coordinates in WGS84 (latitude and longitude)
%
%   Notes:
%
%   This function only checks for the input to be N by 2, not the format, 
%   order, nor consistency of the input coordinates.
%
%   Example - SAFER Lindholmen Göteborg
%
%       x = 6404638;
%       y = 1269670;
%       [lat_and_long] = RT90_2_WGS84([x,y])     
%
% Copyright 2013 Chalmers.
% Marco Dozza - Chalmers - 130318 - marco.dozza@chalmers.se
% The algorithm used in this function is adapted from: https://gist.github.com/moonhouse/2602783

function [coord_WGS84] = RT90_2_WGS84 (varargin)

switch nargin
    case 0
        error('No input given')
    case 1
        if size(varargin{1},2) == 2
            x = varargin{1}(:,1);
            y = varargin{1}(:,2);
        else
            error('Input should be an N by 2 matrix');
        end
    otherwise
        error('To many input arguments: this function expects an N by 2 matrix with RT90 coordinates')
end

Coord_WGS84 = [];

deg_to_rad = pi / 180;
axis = 6378137;
flattening = 1 / 298.257222101 ;
central_meridian = 15 + 48/60 + 22.624306/3600;
scale = 1.00000561024;
false_northing = -667.711;
false_easting = 1500064.274;

e2 = flattening * (2 - flattening);
n = flattening / (2 - flattening);
a_roof = axis / (1 + n) * (1 + n*n/4 + n*n*n*n/64);
delta1 = n/2 - 2*n*n/3 + 37*n*n*n/96 - n*n*n*n/360;
delta2 = n*n/48 + n*n*n/15 - 437*n*n*n*n/1440;
delta3 = 17*n*n*n/480 - 37*n*n*n*n/840;
delta4 = 4397*n*n*n*n/161280;

astar = e2 + e2*e2 + e2*e2*e2 + e2*e2*e2*e2;
bstar = -(7*e2*e2 + 17*e2*e2*e2 + 30*e2*e2*e2*e2) / 6;
cstar = (224*e2*e2*e2 + 889*e2*e2*e2*e2) / 120;
dstar = -(4279*e2*e2*e2*e2) / 1260;

lambda_zero = central_meridian * deg_to_rad;
xi = (x - false_northing) / (scale * a_roof);
eta = (y - false_easting) / (scale * a_roof);
xi_prim = xi -...
    delta1*sin(2*xi) .* cosh(2*eta) -...
    delta2*sin(4*xi) .* cosh(4*eta) -...
    delta3*sin(6*xi) .* cosh(6*eta) -...
    delta4*sin(8*xi) .* cosh(8*eta);
eta_prim = eta -...
    delta1*cos(2*xi) .* sinh(2*eta) -...
    delta2*cos(4*xi) .* sinh(4*eta) -...
    delta3*cos(6*xi) .* sinh(6*eta) -...
    delta4*cos(8*xi) .* sinh(8*eta);

phi_star = asin(sin(xi_prim) ./ cosh(eta_prim));
delta_lambda = atan(sinh(eta_prim) ./ cos(xi_prim));

lon_radian = lambda_zero + delta_lambda;
lat_radian = phi_star + sin(phi_star) .* cos(phi_star) .* (astar +...
    bstar*(sin(phi_star).^2) +...
    cstar*(sin(phi_star).^4) +...
    dstar*(sin(phi_star).^6));

coord_WGS84 = [(lat_radian * 1800000 / pi)/10000,(lon_radian * 1800000 / pi)/10000];

