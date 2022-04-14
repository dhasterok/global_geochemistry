function H = computehp(K2O,U,Th)
% COMPUTEHP - Computes heat production per unit mass
%
%   H = computehp(K2O,U,Th) computes the heat production per unit mass with
%   K2O concentration in wt.% and U and Th in ppm.
%
%   To compute heat production, A = rho*H, where rho is the density
%

% modified from Rybach (1988) for use with K2O instead of K
H = 10^-5*(3.48 * 2 * molecularwt('K')/molecularwt('K2O') * abs(K2O) ...
    + 9.52 * abs(U) ...
    + 2.56 * abs(Th));

return
