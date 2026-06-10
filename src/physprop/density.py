"""
Rock density estimation from whole-rock geochemistry or seismic velocity.

Two independent methods are implemented:

1. Christensen & Mooney (1995): density from P-wave velocity.
2. Empirical oxide model (Hasterok et al.): density from geochemical indices.

References:
    Christensen, N.I. & Mooney, W.D. (1995) Seismic velocity structure and
    composition of the continental crust: a global view. Journal of Geophysical
    Research 100(B7), 9761–9788.

    Behn, M.D. & Kelemen, P.B. (2003) doi:10.1029/2002GC000393
"""

import numpy as np
import pandas as pd


# Shift applied to the geochemical oxide model (calibration offset, kg m⁻³)
_DENSITY_SHIFT = -120.0


def density_cm(vp):
    """Estimate density (kg m⁻³) from P-wave velocity using Christensen & Mooney (1995).

    Parameters
    ----------
    vp : array-like
        P-wave velocity in km/s.

    Returns
    -------
    numpy.ndarray
        Density in kg m⁻³.
    """
    vp = np.asarray(vp, dtype=float)
    return 4929.0 - 13294.0 / vp


def density_geochem(data, shift=_DENSITY_SHIFT, inplace=False):
    """Estimate density (kg m⁻³) from geochemical indices.

    Uses the empirical regression:
        ρ = shift + 2606.5 + 174.7×Fe_number − 12.0×MALI
              + 49.6×ASI + 636.0×maficity

    Requires geochemical indices (added by ``geochem_index()``).
    Outliers outside [2400, 3600] kg m⁻³ are set to NaN.

    Parameters
    ----------
    data : pandas.DataFrame
        Must contain 'Fe_number', 'MALI', 'ASI', 'maficity' columns, plus
        normalized major oxides (for quality filtering).
    shift : float
        Calibration offset in kg m⁻³.  Default −120.
    inplace : bool

    Returns
    -------
    pandas.DataFrame
        With column 'density_bk' (kg m⁻³) added.
    """
    from src.classification.indices import geochem_index
    from src.processing.normalize import oxide_norm, DEFAULT_OXIDES

    if not inplace:
        data = data.copy()

    # Normalize to the 9-oxide system used in the original calibration
    _ox9 = ['sio2', 'tio2', 'al2o3', 'feo_tot', 'mgo', 'cao', 'na2o', 'k2o', 'p2o5']
    tmp = oxide_norm(data, oxides=_ox9, normalization='anhydrous', verbose=False)
    tmp = geochem_index(tmp)

    rho = (shift + 2606.5
           + 174.7 * tmp['Fe_number']
           - 12.0  * tmp['MALI']
           + 49.6  * tmp['ASI']
           + 636.0 * tmp['maficity'])

    # Reject samples where any individual oxide exceeds 100 wt% (bad normalization)
    bad_ox = (tmp[_ox9] > 100).any(axis=1)
    rho[bad_ox] = np.nan

    # Reject physically unreasonable densities
    rho[(rho < 2400) | (rho > 3600)] = np.nan

    data['density_bk'] = rho.values
    return data


def densest(data, inplace=False):
    """Estimate density by both methods and pick the best available estimate.

    Runs ``vpest``, ``density_cm``, and ``density_geochem`` in sequence.
    'density_model' is set to the geochemical model (density_bk) where
    available, falling back to the Christensen–Mooney estimate (density_cm).

    Parameters
    ----------
    data : pandas.DataFrame
    inplace : bool

    Returns
    -------
    pandas.DataFrame
        Adds 'p_velocity', 'density_bk', 'density_cm', 'density_model'.
    """
    from src.physprop.seismic import vpest

    if not inplace:
        data = data.copy()

    data = vpest(data)
    data['density_cm_val'] = density_cm(data['p_velocity'])
    data = density_geochem(data, inplace=True)

    # Best estimate: prefer geochemical model, fall back to Christensen-Mooney
    data['density_model'] = data['density_bk'].where(
        data['density_bk'].notna(), data['density_cm_val']
    )
    return data
