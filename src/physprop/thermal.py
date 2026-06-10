"""
Thermal conductivity estimation from whole-rock geochemistry.

Implements the Willcocks et al. empirical model based on oxide composition.
Valid range: 1.5 – 5 W m⁻¹ K⁻¹.
"""

import numpy as np
import pandas as pd


_TC_OXIDES = ['sio2', 'tio2', 'al2o3', 'feo_tot', 'mgo', 'cao', 'na2o', 'k2o']

# Regression coefficients from Willcocks et al.
_TC_COEFF = {
    'sio2':    1.5542,
    'tio2':   -6.5211,
    'al2o3':  -0.2547,
    'feo_tot': 2.0660,
    'mgo':     0.5472,
    'cao':     0.2258,
    'na2o':   -2.4017,
    'k2o':    -0.6721,
}

# Valid range for thermal conductivity (W m⁻¹ K⁻¹)
_TC_MIN = 1.5
_TC_MAX = 5.0


def tcest(data, inplace=False):
    """Estimate thermal conductivity (W m⁻¹ K⁻¹) from oxide composition.

    Uses the Willcocks et al. model:

        κ = exp(Σᵢ aᵢ × Xᵢ / Σ Xᵢ)

    where Xᵢ are the major oxides normalized anhydrously and aᵢ are
    fitted coefficients.  Values outside [1.5, 5] W m⁻¹ K⁻¹ are set to NaN.

    Parameters
    ----------
    data : pandas.DataFrame
        Major oxide columns in wt%.
    inplace : bool

    Returns
    -------
    pandas.DataFrame
        With column 'thermal_conductivity' (W m⁻¹ K⁻¹) added.
    """
    from src.processing.normalize import oxide_norm

    if not inplace:
        data = data.copy()

    tmp = oxide_norm(data, oxides=_TC_OXIDES, normalization='anhydrous',
                     verbose=False)

    ox_vals = tmp[_TC_OXIDES].fillna(0.0)
    total = ox_vals.sum(axis=1).replace(0, np.nan)

    numerator = sum(
        _TC_COEFF[ox] * ox_vals[ox] for ox in _TC_OXIDES
    )

    tc = np.exp(numerator / total)
    tc[(tc < _TC_MIN) | (tc > _TC_MAX)] = np.nan

    data['thermal_conductivity'] = tc.values
    return data
