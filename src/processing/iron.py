"""
Iron oxidation state consolidation for whole-rock geochemical data.

Databases commonly report iron in multiple forms: FeO (ferrous), Fe2O3
(ferric), Fe2O3_tot (total iron as ferric), or FeO_tot (total iron as ferrous).
This module converts whichever combination is present to a single 'feo_tot'
column and, where both FeO and Fe2O3 are measured, records the Fe2+/ΣFe ratio.

Conversion factor:  FeO_tot = Fe2O3 × 2 × MW(FeO) / MW(Fe2O3)
                             = Fe2O3 × 0.8998...

Reference:
    Middlemost (1989) Iron oxidation ratios, norms and the classification of
    volcanic rocks.  Chemical Geology 77, 19–26.
"""

import numpy as np
import pandas as pd
from src.utils.molecular import MolecularWeightCalculator

_mwc = MolecularWeightCalculator()

# Conversion: Fe2O3 (wt%) → FeO (wt%)  [multiply Fe2O3 by this]
_FE_FACTOR = 2 * _mwc.molecular_weight('FeO') / _mwc.molecular_weight('Fe2O3')


def fefix(data, inplace=False):
    """Consolidate iron reporting to a single 'feo_tot' column.

    Handles six possible iron reporting scenarios in the input data:

    1. Only FeO present          → feo_tot = feo; drop feo
    2. Only Fe2O3 present        → feo_tot = fe2o3 × factor; drop fe2o3
    3. Only Fe2O3_tot present     → feo_tot = fe2o3_tot × factor; drop fe2o3_tot
    4. Both FeO and Fe2O3        → feo_tot = feo + fe2o3 × factor;
                                   record fe2_fe_tot ratio; drop both
    5. Fe2O3_tot, no FeO/feo_tot → feo_tot = fe2o3_tot × factor; drop fe2o3_tot
    6. feo_tot already present   → no conversion, just clean up

    Also handles 'feto3' (USGS synonym for Fe2O3).

    Parameters
    ----------
    data : pandas.DataFrame
        Table with lowercase iron column names.
    inplace : bool
        Modify *data* in place.  Default False (returns a copy).

    Returns
    -------
    pandas.DataFrame
        DataFrame with 'feo_tot' populated and intermediate columns removed.
        Adds 'fe2_fe_tot' (Fe²⁺/ΣFe molar ratio) where both FeO and Fe2O3
        were measured.

    Notes
    -----
    To recover FeO and Fe2O3 from the output::

        import src.utils.molecular as mol
        data['feo']   = mol.mw('FeO')  * data['fe2_fe_tot'] * data['feo_tot']
        data['fe2o3'] = 0.5 * mol.mw('Fe2O3') * (1 - data['fe2_fe_tot']) * data['feo_tot']
    """
    if not inplace:
        data = data.copy()

    cols = set(data.columns)

    # USGS datasets sometimes use 'feto3' as the Fe2O3 column name
    if 'feto3' in cols:
        data['fe2o3'] = data['feto3']
        data.drop(columns=['feto3'], inplace=True)
        cols = set(data.columns)

    has_feo     = 'feo'      in cols
    has_fe2o3   = 'fe2o3'   in cols
    has_fe2o3t  = 'fe2o3_tot' in cols
    has_feot    = 'feo_tot' in cols

    # --- Fast exit: no iron columns at all ---
    if not any([has_feo, has_fe2o3, has_fe2o3t, has_feot]):
        return data

    # --- Case 1: only feo ---
    if has_feo and not has_fe2o3 and not has_fe2o3t and not has_feot:
        data.rename(columns={'feo': 'feo_tot'}, inplace=True)
        return data

    # --- Case 2: only fe2o3 ---
    if has_fe2o3 and not has_feo and not has_feot and not has_fe2o3t:
        data['feo_tot'] = data['fe2o3'] * _FE_FACTOR
        data.drop(columns=['fe2o3'], inplace=True)
        return data

    # --- Case 3: only fe2o3_tot ---
    if has_fe2o3t and not has_feo and not has_feot and not has_fe2o3:
        data['feo_tot'] = data['fe2o3_tot'] * _FE_FACTOR
        data.drop(columns=['fe2o3_tot'], inplace=True)
        return data

    # --- General case: one or more columns coexist ---
    if not has_feot:
        data['feo_tot'] = np.nan
    if not has_fe2o3t:
        data['fe2o3_tot'] = np.nan

    data['fe2_fe_tot'] = np.nan

    # Case 4a: feo only, no fe2o3 → direct assign
    mask = (
        data['feo'].gt(0) &
        data['feo_tot'].isna() &
        (~data['fe2o3'].gt(0) if has_fe2o3 else True)
    ) if has_feo else pd.Series(False, index=data.index)
    data.loc[mask, 'feo_tot'] = data.loc[mask, 'feo']

    # Case 4b: both feo and fe2o3 → combine and record ratio
    if has_feo and has_fe2o3:
        mask = data['fe2o3'].gt(0) & data['feo'].gt(0)
        feo_val   = data.loc[mask, 'feo']
        fe2o3_val = data.loc[mask, 'fe2o3']
        data.loc[mask, 'feo_tot'] = feo_val + fe2o3_val * _FE_FACTOR

        # Molar Fe2+/ΣFe ratio
        mol_feo   = feo_val   / _mwc.molecular_weight('FeO')
        mol_fe2o3 = fe2o3_val / _mwc.molecular_weight('Fe2O3')
        data.loc[mask, 'fe2_fe_tot'] = mol_feo / (mol_feo + 2 * mol_fe2o3)

    # Case 5: fe2o3_tot where fe2o3 alone present with no feo
    if has_fe2o3:
        mask = (
            data['fe2o3'].gt(0) &
            data['fe2o3_tot'].isna() &
            (~data['feo'].gt(0) if has_feo else True)
        )
        data.loc[mask, 'fe2o3_tot'] = data.loc[mask, 'fe2o3']

    # Case 6: fe2o3_tot → feo_tot where feo_tot still missing
    mask = data['fe2o3_tot'].notna() & data['feo_tot'].isna()
    data.loc[mask, 'feo_tot'] = data.loc[mask, 'fe2o3_tot'] * _FE_FACTOR

    # Clean up intermediate columns
    for col in ['feo', 'fe2o3']:
        if col in data.columns:
            data.drop(columns=[col], inplace=True)

    return data


def fe_ratio(feo, fe2o3):
    """Compute Fe2+/ΣFe molar ratio from FeO and Fe2O3 wt%.

    Parameters
    ----------
    feo, fe2o3 : array-like
        FeO and Fe2O3 concentrations in wt%.

    Returns
    -------
    numpy.ndarray
        Fe2+/(Fe2+ + Fe3+) molar ratio.
    """
    feo   = np.asarray(feo,   dtype=float)
    fe2o3 = np.asarray(fe2o3, dtype=float)
    mol_feo   = feo   / _mwc.molecular_weight('FeO')
    mol_fe2o3 = fe2o3 / _mwc.molecular_weight('Fe2O3')
    with np.errstate(invalid='ignore', divide='ignore'):
        ratio = mol_feo / (mol_feo + 2 * mol_fe2o3)
    ratio[~np.isfinite(ratio)] = np.nan
    return ratio
