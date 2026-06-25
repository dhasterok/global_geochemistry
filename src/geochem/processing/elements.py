"""
Element summation and grouping utilities for whole-rock geochemical data.

Element concentrations in trace element columns are assumed to be in ppm
(column names end in '_ppm').  Sums follow the conventions of Taylor &
McLennan (1985) and Sun & McDonough (1989) for REE groupings.
"""

import numpy as np
import pandas as pd


# -------------------------------------------------------------------------
# REE definitions
# -------------------------------------------------------------------------

_LREE_FULL    = ['la_ppm', 'ce_ppm', 'pr_ppm', 'nd_ppm', 'sm_ppm', 'eu_ppm', 'gd_ppm']
_HREE_FULL    = ['tb_ppm', 'dy_ppm', 'ho_ppm', 'er_ppm', 'tm_ppm', 'yb_ppm', 'lu_ppm']
_REE_FULL     = _LREE_FULL + _HREE_FULL

_LREE_REDUCED = ['la_ppm', 'ce_ppm', 'nd_ppm']
_HREE_REDUCED = ['er_ppm', 'yb_ppm', 'lu_ppm']
_MREE         = ['sm_ppm', 'eu_ppm', 'gd_ppm']
_REE_REDUCED  = _LREE_REDUCED + _MREE + _HREE_REDUCED


def sumree(data, scheme='full', inplace=False):
    """Sum rare earth element concentrations into group totals.

    Sums only detected values (> 0); below-detection-limit (BDL) and NaN
    entries are excluded from the running sum.  If *all* elements in a
    group are BDL or absent, the group sum is set to NaN.

    Parameters
    ----------
    data : pandas.DataFrame
        Must contain ppm columns named e.g. ``'la_ppm'``, ``'ce_ppm'``, …
    scheme : {'full', 'reduced'}
        'full' (default): uses all 14 REE (La–Lu).
        'reduced': uses a subset (La, Ce, Nd | Sm, Eu, Gd | Er, Yb, Lu) to
        reduce the impact of BDL values.
    inplace : bool
        Modify *data* in place.  Default False.

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns ``'ree'``, ``'lree'``, ``'mree'``, ``'hree'``
        added (total, light, middle, heavy REE sums in ppm).

    Notes
    -----
    The MATLAB original had a logic bug in the BDL-masking step
    (``all(X > 0)`` operated column-wise instead of row-wise).  This
    implementation uses the correct row-wise check.

    References
    ----------
    Taylor & McLennan (1985) The continental crust: its composition and
    evolution. Blackwell.
    """
    if scheme == 'full':
        lree_cols = _LREE_FULL
        hree_cols = _HREE_FULL
        ree_cols  = _REE_FULL
        print('Computing REE sums using full element set.')
    elif scheme == 'reduced':
        lree_cols = _LREE_REDUCED
        hree_cols = _HREE_REDUCED
        ree_cols  = _REE_REDUCED
        print('Computing REE sums using reduced element set.')
    else:
        raise ValueError(f"scheme must be 'full' or 'reduced', got '{scheme}'.")

    if not inplace:
        data = data.copy()

    def _group_sum(df, cols):
        present = [c for c in cols if c in df.columns]
        if not present:
            return pd.Series(np.nan, index=df.index)
        vals = df[present].copy()
        vals[vals <= 0] = np.nan          # treat BDL (≤ 0) as missing
        s = vals.sum(axis=1, min_count=1) # NaN if all missing
        return s

    data['ree']  = _group_sum(data, ree_cols)
    data['lree'] = _group_sum(data, lree_cols)
    data['mree'] = _group_sum(data, _MREE)
    data['hree'] = _group_sum(data, hree_cols)

    return data


# -------------------------------------------------------------------------
# K2O utilities
# -------------------------------------------------------------------------

def get_k2o(data):
    """Return K2O column, preferring measured over converted values.

    Looks for 'k2o' first, then 'k_ppm' converted to oxide.  Returns a
    Series with NaN where neither is available.

    Parameters
    ----------
    data : pandas.DataFrame

    Returns
    -------
    pandas.Series
        K2O in wt%.
    """
    from global_geochemistry.utils.molecular import MolecularWeightCalculator
    _mwc = MolecularWeightCalculator()

    if 'k2o' in data.columns:
        return data['k2o'].copy()
    if 'k_ppm' in data.columns:
        # ppm K → wt% K2O:  K_ppm × MW(K2O) / (2 × MW(K) × 1e4)
        factor = _mwc.molecular_weight('K2O') / (2 * _mwc.molecular_weight('K') * 1e4)
        return (data['k_ppm'] * factor).where(data['k_ppm'] > 0)
    return pd.Series(np.nan, index=data.index)


def fix_k(data, inplace=False):
    """Fill missing K2O from K_ppm where K2O is absent or non-positive.

    Parameters
    ----------
    data : pandas.DataFrame
    inplace : bool

    Returns
    -------
    pandas.DataFrame
    """
    if not inplace:
        data = data.copy()
    if 'k_ppm' not in data.columns:
        return data
    if 'k2o' not in data.columns:
        data['k2o'] = np.nan

    from global_geochemistry.utils.molecular import MolecularWeightCalculator
    _mwc = MolecularWeightCalculator()
    factor = _mwc.molecular_weight('K2O') / (2 * _mwc.molecular_weight('K') * 1e4)
    k2o_from_ppm = data['k_ppm'] * factor

    fill_mask = (data['k2o'].isna() | (data['k2o'] <= 0)) & (data['k_ppm'] > 0)
    data.loc[fill_mask, 'k2o'] = k2o_from_ppm[fill_mask]
    return data
