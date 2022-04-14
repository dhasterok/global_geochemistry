function data = vpest(data)
% estimation of seismic velocity from Behn & Kelemen [G3,
% 2003, doi:10.1029/2002GC000393]

% Behn & Kelemen use 7 oxides.  Normalize compositions to these 7 for
% density and velocity estimates.
oxides = {'SiO2'; 'Al2O3'; 'FeO'; 'MgO'; 'CaO'; 'Na2O'; 'K2O'};
tmp = oxide_norm(data,oxides);

Vp = 6.9 ...
- 0.011 * tmp.sio2 ...
+ 0.037 * tmp.mgo ...
+ 0.045 * tmp.cao;

data.p_velocity = Vp;

return
