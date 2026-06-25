"""
Seismic P-wave velocity estimation from whole-rock geochemistry.

Implements the Behn & Kelemen (2003) empirical regression of Vp against
major oxide composition, calibrated for lower-crustal and mantle lithologies.

Reference:
    Behn, M.D. & Kelemen, P.B. (2003) Relationship between seismic P-wave
    velocity and the composition of anhydrous igneous and meta-igneous rocks.
    Geochemistry, Geophysics, Geosystems 4(5), 1041.
    doi:10.1029/2002GC000393
"""

import numpy as np
import pandas as pd


# 7-oxide system used by Behn & Kelemen (2003)
_BK_OXIDES = ['sio2', 'al2o3', 'feo_tot', 'mgo', 'cao', 'na2o', 'k2o']


def vpest(data, inplace=False):
    """Estimate P-wave velocity (km/s) from oxide composition.

    Uses the Behn & Kelemen (2003) regression on oxides normalized to
    the 7-component SiO2–Al2O3–FeO–MgO–CaO–Na2O–K2O system:

        Vp = 6.9 − 0.011×SiO2 + 0.037×MgO + 0.045×CaO

    Parameters
    ----------
    data : pandas.DataFrame
        Must contain major oxide columns (wt%, normalized).
    inplace : bool

    Returns
    -------
    pandas.DataFrame
        With column 'p_velocity' (km/s) added.
    """
    from global_geochemistry.geochem.processing.normalize import oxide_norm

    if not inplace:
        data = data.copy()

    tmp = oxide_norm(data, oxides=_BK_OXIDES, normalization='anhydrous',
                     verbose=False)

    sio2 = tmp['sio2'].fillna(0)
    mgo  = tmp['mgo'].fillna(0)
    cao  = tmp['cao'].fillna(0)

    vp = 6.9 - 0.011 * sio2 + 0.037 * mgo + 0.045 * cao

    # Mask samples where normalization failed (all oxides absent)
    no_data = tmp[_BK_OXIDES].isna().all(axis=1)
    vp[no_data] = np.nan

    data['p_velocity'] = vp.values
    return data


def vsest(vp):
    """Estimate S-wave velocity (km/s) from P-wave velocity.

    Uses the empirical Vp–Vs relationship for crustal rocks:
        Vs ≈ Vp / 1.73   (Poisson's ratio ≈ 0.25)

    Parameters
    ----------
    vp : array-like
        P-wave velocity in km/s.

    Returns
    -------
    numpy.ndarray
    """
    return np.asarray(vp, dtype=float) / 1.73
