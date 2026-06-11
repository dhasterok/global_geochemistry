function data = vsest(data)
% estimation of seismic velocity from Behn & Kelemen [G3,
% 2003, doi:10.1029/2002GC000393]... fitting done by Matt Linke (Honours
% 2020)

% Behn & Kelemen use 7 oxides.  Normalize compositions to these 7 for
% density and velocity estimates.
oxides = {'SiO2'; 'Al2O3'; 'FeO'; 'MgO'; 'CaO'; 'Na2O'; 'K2O'};
tmp = oxide_norm(data,'Normalization','anhydrous','Oxides',oxides);

Vs = 4.5959 ...
    - 0.0091 * tmp.sio2 ...
    - 0.0150 * tmp.al2o3 ...
    + 0.0142 * tmp.mgo;

data.s_velocity = Vs;

return