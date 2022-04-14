function data = vpest2(data)
% estimation of seismic velocity from Behn & Kelemen [G3,
% 2003, doi:10.1029/2002GC000393]
%Without normalisation

Vp = 6.9 ...
    - 0.011 * tmp.sio2 ...
    + 0.037 * tmp.mgo ...
    + 0.045 * tmp.cao;

data.p_velocity = Vp;

return
