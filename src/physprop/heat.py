"""
Radiogenic heat production estimation from whole-rock geochemistry.

Heat production is calculated from K, U, and Th concentrations using the
Rybach (1988) formulation.  Volumetric heat production requires a density
estimate (see src/physprop/density.py).

References:
    Rybach, L. (1988) Determination of heat production rate. In: Haenel R.,
    Rybach L., Stegena L. (eds) Handbook of Terrestrial Heat-Flow Density
    Determination. Kluwer, Dordrecht. pp. 125–142.
"""

import numpy as np
import pandas as pd
from src.utils.molecular import MolecularWeightCalculator

_mwc = MolecularWeightCalculator()

# Rybach (1988) constants (μW kg⁻¹ per unit concentration)
# K: per wt% K  (note: MATLAB uses K2O; we convert internally)
# U: per ppm U
# Th: per ppm Th
_CK  = 3.48     # μW kg⁻¹ per wt% K
_CU  = 9.52     # μW kg⁻¹ per ppm U
_CTH = 2.56     # μW kg⁻¹ per ppm Th
_SCALE = 1e-5   # overall scale factor (result in μW kg⁻¹)

# Conversion from K2O wt% to K wt%
_K2O_TO_K = 2 * _mwc.molecular_weight('K') / _mwc.molecular_weight('K2O')


def compute_hp(k2o, u, th):
    """Compute specific heat production (μW kg⁻¹) from K2O, U, and Th.

    Parameters
    ----------
    k2o : array-like
        K2O concentration in wt%.
    u : array-like
        Uranium concentration in ppm.
    th : array-like
        Thorium concentration in ppm.

    Returns
    -------
    numpy.ndarray
        Heat production in μW kg⁻¹.

    Notes
    -----
    Formula (Rybach 1988):
        H = 1e-5 × (3.48 × [K]_wt% + 9.52 × [U]_ppm + 2.56 × [Th]_ppm)

    where [K]_wt% = 2 × MW(K)/MW(K2O) × [K2O]_wt%.
    """
    k2o = np.asarray(k2o, dtype=float)
    u   = np.asarray(u,   dtype=float)
    th  = np.asarray(th,  dtype=float)
    k   = _K2O_TO_K * np.abs(k2o)
    return _SCALE * (_CK * k + _CU * np.abs(u) + _CTH * np.abs(th))


def hpest(data, bdl_k2o=-0.01, bdl_u=-1.0, bdl_th=-1.0, verbose=True):
    """Estimate heat production and add it to a geochemical data table.

    Adds columns 'heat_production_mass' (μW kg⁻¹) and 'heat_production'
    (μW m⁻³, if density is available).  Negative values signal that at
    least one input was below the detection limit.

    Parameters
    ----------
    data : pandas.DataFrame
        Must contain 'u_ppm', 'th_ppm', and either 'k2o' or 'k_ppm'.
    bdl_k2o, bdl_u, bdl_th : float
        Below-detection-limit substitute values.  Zeros in the source data
        are replaced with these values before computing; negative substitutes
        indicate BDL.  Defaults: K2O = -0.01 wt%, U = Th = -1 ppm.
    verbose : bool

    Returns
    -------
    pandas.DataFrame
        Modified copy of *data*.
    """
    data = data.copy()

    # Resolve K2O from k2o or k_ppm
    if 'k2o' in data.columns:
        k2o = data['k2o'].copy()
    elif 'k_ppm' in data.columns:
        factor = _mwc.molecular_weight('K2O') / (2 * _mwc.molecular_weight('K') * 1e4)
        k2o = data['k_ppm'] * factor
    else:
        k2o = pd.Series(np.nan, index=data.index)

    u  = data['u_ppm'].copy()  if 'u_ppm'  in data.columns else pd.Series(np.nan, index=data.index)
    th = data['th_ppm'].copy() if 'th_ppm' in data.columns else pd.Series(np.nan, index=data.index)

    if verbose:
        print(f'BDL substitutes:  K2O = {bdl_k2o} wt%,  U = {bdl_u} ppm,  Th = {bdl_th} ppm')

    # Replace zero concentrations with BDL substitute
    k2o = k2o.where(k2o != 0, bdl_k2o)
    u   = u.where(u   != 0, bdl_u)
    th  = th.where(th  != 0, bdl_th)

    # Sign flag: +1 all detected, -1 any BDL
    all_detected = (k2o > 0) & (u > 0) & (th > 0)
    sign = np.where(all_detected, 1.0, -1.0)
    sign = pd.Series(sign, index=data.index, dtype=float)
    sign[k2o.isna() | u.isna() | th.isna()] = np.nan

    data['heat_production_mass'] = sign * compute_hp(
        k2o.abs(), u.abs(), th.abs()
    )

    # Volumetric heat production (μW m⁻³)
    if 'density_model' not in data.columns:
        from src.physprop.seismic import vpest
        data = vpest(data)
        from src.physprop.density import density_cm
        data['density_model'] = density_cm(data['p_velocity'])

    data['heat_production'] = data['heat_production_mass'] * data['density_model']

    return data
