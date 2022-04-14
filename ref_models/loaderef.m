function eref = loaderef;
% LOADEREF - Loads reference geochemistry reservoirs
%
%   eref = loaderef will load reference geochemistry reservoirs into the
%   table eref.  The physical properties will also be computed for the
%   reference reservoirs.

eref = readtable('earthref.xlsx');
%eref = ecread2('earthref.xlsx');

eref = fefix(eref);

% compute geochemical indices
eref = geochem_index(eref);

disp('Computing rock properties...');
% estimate seismic velocities
eref = vpest(eref);

% estimate density
eref = densest2(eref);

% estimate heat production
eref = hpest(eref);

% estimate thermal conductivity
eref = tcest(eref);

eref.sigma = logical(eref.sigma);

eref.p_velocity(eref.sigma) = NaN;
eref.s_velocity(eref.sigma) = NaN;
eref.heat_production_mass(eref.sigma) = NaN;
eref.heat_production(eref.sigma) = NaN;
eref.density_bk(eref.sigma) = NaN;

return
