function data = densest2(data)
% DENSEST - estimates density using several methods
%
% data = densest2(data)
% Density is estimated using two methods:
%    (1) Christensen and Mooney (1995) where density is estimated from sesimic
%    velocity; and
%    (2) a model using a empirical density equation as a function of oxide
%    content.

% Behn & Kelemen use 7 oxides.  Normalize compositions to these 7 for
% density and velocity estimates.
oxides = {'SiO2'; 'TiO2'; 'Al2O3'; 'FeO'; 'MgO'; 'CaO'; 'Na2O'; 'K2O'; 'P2O5'};
tmp = oxide_norm(data,oxides);
tmp = geochem_index(tmp);

% Indicies for various density models
% ----------------------------------------------
igind = rockgroup(data,'igneous protolith');
sedind = rockgroup(data,'sedimentary protolith');

% Group 3 - Igneous carbonatites
carbonatite = strcmp('carbonatite',data.rock_type);
% Group 4 - Sedimentary carbonates
carbonate = strcmp('dolomite',data.rock_type) | strcmp('limestone',data.rock_type);
% Group 1 - Low-magnesian, igneous and sedimentary rocks
%   sans carbonates/carbonatites
low_mg = ((data.mgo < 18 & rockgroup(data,'all igneous')) | rockgroup(data,'all seds')) & ~(carbonatite | carbonate);
% Group 2 - High-magnesian, igneous and sedimentary rocks
%   sans carbonates/carbonatites
high_mg = (data.mgo >= 18 & rockgroup(data,'all igneous')) & ~(carbonatite | carbonate);


% Estimate density from seismic velocity
% ----------------------------------------------
% estimation of density using Christensen and Mooney [JGR, 1995]
data.density_cm = 4929 - 13294./data.p_velocity;


% Estimate density from thermodynamic calculations
% ----------------------------------------------
% Note Behn and Kelemen (2003) did not produce a model in their paper, but did
% compute density as part of their velocity caluclations.

shift = -120;
rho = shift + 2606.5 + 174.7*tmp.Fe_number - 12.0*tmp.MALI + 49.6*tmp.ASI + 636.0*tmp.maficity;

% only use data with reasonable oxide percentages
ind = (tmp.sio2 > 100 | ...
    tmp.tio2 > 100 | ...
    tmp.al2o3 > 100 | ...
    tmp.feo_tot > 100 | ...
    tmp.mgo > 100 | ...
    tmp.cao > 100 | ...
    tmp.na2o > 100 | ...
    tmp.k2o > 100 | ...
    tmp.p2o5 > 100);

rho(ind) = NaN;

% remove outliers
rho(rho<2400 | rho>3600) = NaN;
rho = rho(:);

data.density_bk = rho;

% Estimate density from empirical fits to data
% there are four different models
% ----------------------------------------------
data.density_sio2 = 2531.463956 + 216.172176*data.Fe_number + ...
    608.428719*data.maficity - 9.948247*data.MALI;
data.density_himg = 3143 - 10.6*data.mgo + 6.2*data.cao;
data.density_igcarb = 2733 + 13.6*data.feo_tot + 8.8*data.p2o5;
data.density_sedcarb = 3402 - 8.0*data.sio2 - 4.5*data.mgo - 7.6*data.cao;

% default model
data.density_model = data.density_bk;

% use a better model if it is possible to compute it
data.density_model(low_mg) = data.density_sio2(low_mg);
data.density_model(high_mg) = data.density_sio2(high_mg);
data.density_model(carbonatite) = data.density_igcarb(carbonatite);
data.density_model(carbonate) = data.density_sedcarb(carbonate);

% restrict density between 2400 and 3600...remember these are zero-porosity
% models so they should be denser than typically estimated for sediments.
data.density_model(data.density_model < 2400 | data.density_model > 3600) = NaN;

return

