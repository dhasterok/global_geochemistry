"""
Oxide normalization and cation conversion for whole-rock geochemical data.

All compositions are in wt% unless noted (ppm fields end in _ppm).
Iron is assumed to be reported as total iron as FeO in the column 'feo_tot'.

References:
    Le Maitre et al. (2002) A classification of igneous rocks and a glossary
    of terms, Cambridge University Press.
"""

import re
import warnings
import numpy as np
import pandas as pd
from src.utils.molecular import MolecularWeightCalculator

_mwc = MolecularWeightCalculator()


def mw(formula):
    """Return molecular weight of a chemical formula (g/mol)."""
    return _mwc.molecular_weight(formula)


# Standard anhydrous oxide list (total iron as FeO)
DEFAULT_OXIDES = [
    'sio2', 'tio2', 'al2o3', 'feo_tot',
    'mgo', 'cao', 'na2o', 'k2o', 'p2o5',
    'cr2o3', 'mno', 'nio',
]

# Volatile species tracked alongside major oxides
VOLATILES = ['h2o_tot', 'co2', 'so3', 'f_ppm', 'cl_ppm']

# Conversion factors: volatile columns reported in ppm need → wt%
_VOLATILE_CF = {
    'h2o_tot': 1.0,
    'co2':     1.0,
    'so3':     1.0,
    'f_ppm':   1e-4,
    'cl_ppm':  1e-4,
}


def oxide_norm(data, oxides=None, normalization='anhydrous', total_tol=None,
               verbose=True):
    """Normalize major oxide concentrations to 100 wt%.

    Parameters
    ----------
    data : pandas.DataFrame
        Table with lowercase oxide column names.  At minimum should contain
        some subset of ``DEFAULT_OXIDES``.
    oxides : list of str, optional
        Oxide columns to include in the normalization sum.  Defaults to
        ``DEFAULT_OXIDES``.  'feo' is silently aliased to 'feo_tot'.
    normalization : {'anhydrous', 'hydrous'}
        'anhydrous' (default) normalizes to the anhydrous oxide sum.
        'hydrous' normalizes to oxide sum + LOI/volatiles.
    total_tol : float or None
        If given, remove samples whose unnormalized oxide total deviates
        more than *total_tol* wt% from 100 (samples without any major oxides
        are kept regardless).  E.g. ``total_tol=10`` keeps only samples with
        totals between 90 and 110 wt%.  Default ``None`` keeps all samples.
    verbose : bool
        Print summary statistics.  Default True.

    Returns
    -------
    pandas.DataFrame
        Copy of *data* with oxide columns scaled to sum to 100, plus added
        columns ``total_ox`` and ``loi``.
    """
    if oxides is None:
        oxides = list(DEFAULT_OXIDES)
    else:
        oxides = [o.lower() for o in oxides]

    # Alias feo → feo_tot
    oxides = ['feo_tot' if o == 'feo' else o for o in oxides]

    # Drop volatiles from oxide list in case caller included them
    oxides = [o for o in oxides if o not in VOLATILES]

    data = data.copy()

    # ------------------------------------------------------------------
    # Ensure all required columns exist
    # ------------------------------------------------------------------
    for col in oxides + ['loi'] + VOLATILES + ['h2o_plus', 'h2o_minus']:
        if col not in data.columns:
            data[col] = np.nan

    # ------------------------------------------------------------------
    # Carbonate conversions (fill CaO/MgO/CO2 where absent)
    # ------------------------------------------------------------------
    if 'caco3' in data.columns:
        mask = (data['caco3'] > 0) & (data['co2'] < 0) & (data['cao'] < 0)
        cao_frac = mw('CaO') / mw('CaCO3')
        co2_frac = mw('CO2') / mw('CaCO3')
        data.loc[mask, 'cao']  = data.loc[mask, 'caco3'] * cao_frac
        data.loc[mask, 'co2']  = data.loc[mask, 'caco3'] * co2_frac
        data.drop(columns=['caco3'], inplace=True)

    if 'mgco3' in data.columns:
        mask = (data['mgco3'] > 0) & (data['co2'] < 0) & (data['mgo'] < 0)
        mgo_frac = mw('MgO') / mw('MgCO3')
        co2_frac = mw('CO2') / mw('MgCO3')
        data.loc[mask, 'mgo']  = data.loc[mask, 'mgco3'] * mgo_frac
        data.loc[mask, 'co2']  = data.loc[mask, 'mgco3'] * co2_frac
        data.drop(columns=['mgco3'], inplace=True)

    # ------------------------------------------------------------------
    # Combine H2O+ and H2O- into H2O_tot where H2O_tot is absent
    # ------------------------------------------------------------------
    missing_h2o = data['h2o_tot'].isna() | (data['h2o_tot'] < 0)
    data.loc[missing_h2o, 'h2o_tot'] = (
        data.loc[missing_h2o, ['h2o_plus', 'h2o_minus']]
        .clip(lower=0).sum(axis=1)
        .replace(0, np.nan)
    )

    # ------------------------------------------------------------------
    # Compute LOI from individual volatiles where LOI is absent
    # ------------------------------------------------------------------
    vol_wt = pd.DataFrame({
        v: data[v].clip(lower=0) * _VOLATILE_CF[v] for v in VOLATILES
    })
    computed_loi = vol_wt.sum(axis=1)
    computed_loi = computed_loi.where(computed_loi > 0, np.nan)

    # If CO2 > 1.2 × LOI, assume the LOI was measured without carbonate
    co2_exceeds = 1.2 * data['loi'].fillna(0) < data['co2'].fillna(0)
    data.loc[co2_exceeds, 'loi'] = data.loc[co2_exceeds, 'loi'].fillna(0) + \
                                    data.loc[co2_exceeds, 'co2'].fillna(0)

    fill_loi = data['loi'].isna() & computed_loi.notna()
    data.loc[fill_loi, 'loi'] = computed_loi[fill_loi]
    data['total_loi'] = fill_loi.astype(float)  # 1 where LOI was computed

    # ------------------------------------------------------------------
    # Compute oxide total (treating negatives as below-detection → 0)
    # ------------------------------------------------------------------
    ox_vals = data[oxides].clip(lower=0)
    data['total_ox'] = ox_vals.sum(axis=1, min_count=1)  # NaN if all NaN

    # ------------------------------------------------------------------
    # Tolerance filter
    # ------------------------------------------------------------------
    n_before = len(data)
    if total_tol is not None:
        anhydrous_ok = (
            (100 - total_tol < data['total_ox']) &
            (data['total_ox'] < 100 + total_tol)
        )
        hydrous_ok = (
            (100 - total_tol < data['total_ox'] + data['loi'].fillna(0)) &
            (data['total_ox'] + data['loi'].fillna(0) < 100 + total_tol)
        )
        no_data = data['total_ox'].isna()
        data = data[anhydrous_ok | hydrous_ok | no_data].copy()
        if verbose:
            print(f'Tolerance ±{total_tol} wt% filter: removed '
                  f'{n_before - len(data)}/{n_before} samples.')

    # ------------------------------------------------------------------
    # Normalization factor
    # ------------------------------------------------------------------
    nf = np.ones(len(data))
    has_ox = data['total_ox'] > 0

    if normalization == 'anhydrous':
        nf = np.where(has_ox, 100.0 / data['total_ox'], 1.0)
    elif normalization == 'hydrous':
        total_hydrous = data['total_ox'] + data['loi'].fillna(0)
        has_h = total_hydrous > 0
        nf = np.where(has_h, 100.0 / total_hydrous, 1.0)
    else:
        raise ValueError(f"normalization must be 'anhydrous' or 'hydrous', "
                         f"got '{normalization}'.")

    # ------------------------------------------------------------------
    # Apply normalization
    # ------------------------------------------------------------------
    data[oxides]    = data[oxides].multiply(nf, axis=0)
    data['total_ox'] = data['total_ox'] * nf
    data['loi']      = data['loi'] * nf
    for v in VOLATILES:
        data[v] = data[v] * nf

    if verbose:
        tx = data['total_ox']
        nonzero = tx[tx != 0]
        print(f'Normalization ({normalization}) complete.  '
              f'Samples: {len(data)}.  '
              f'Oxide total <40: {(nonzero < 40).sum()}, '
              f'<90: {(nonzero < 90).sum()}, '
              f'>110: {(nonzero > 110).sum()}.')

    return data


def cation_to_oxide(data, inplace=False):
    """Fill missing oxide columns using element concentrations in ppm.

    Where a ppm element column (e.g. ``'si_ppm'``) exists but the
    corresponding oxide column (``'sio2'``) is missing or NaN, convert
    and fill.  Only fills where the oxide is absent — existing oxide data
    are not overwritten.

    Parameters
    ----------
    data : pandas.DataFrame
    inplace : bool
        Modify *data* in place.  Default False (returns a copy).

    Returns
    -------
    pandas.DataFrame
    """
    # Element → oxide mapping: (element_symbol, oxide_formula)
    # Order matters: Fe goes to feo_tot
    _CONVERSIONS = [
        ('Si', 'SiO2'),
        ('Ti', 'TiO2'),
        ('Al', 'Al2O3'),
        ('Cr', 'Cr2O3'),
        ('Fe', 'FeO'),       # → stored in feo_tot
        ('Mg', 'MgO'),
        ('Ca', 'CaO'),
        ('Mn', 'MnO'),
        ('Ni', 'NiO'),
        ('Ba', 'BaO'),
        ('Sr', 'SrO'),
        ('Na', 'Na2O'),
        ('K',  'K2O'),
        ('P',  'P2O5'),
    ]

    if not inplace:
        data = data.copy()

    for el_sym, oxide_formula in _CONVERSIONS:
        el_col  = el_sym.lower() + '_ppm'
        ox_col  = oxide_formula.lower()

        # FeO is stored as feo_tot
        if ox_col == 'feo':
            ox_col = 'feo_tot'

        if el_col not in data.columns:
            continue

        if ox_col not in data.columns:
            data[ox_col] = np.nan

        # Number of cations per formula unit (e.g. Al2O3 → 2 Al)
        n_cat = _count_cations(el_sym, oxide_formula)
        factor = mw(oxide_formula) / (1e4 * n_cat * mw(el_sym))

        converted = data[el_col] * factor
        converted = converted.where(converted <= 100, np.nan)  # reject > 100 wt%

        # Only fill where oxide is absent or negative
        fill_mask = converted > 0
        if (data[ox_col] > 0).any():
            fill_mask = fill_mask & ~(data[ox_col] > 0)

        data.loc[fill_mask, ox_col] = converted[fill_mask]

    return data


def cation_mole_norm(data, oxides, basis=100):
    """Normalize oxide compositions to moles of cations.

    Converts each oxide wt% to moles of the corresponding cation, then
    scales so the cation mole sum equals *basis*.  Adds cation columns
    (element symbol, uppercase first letter) to the returned DataFrame.

    Parameters
    ----------
    data : pandas.DataFrame
        Must contain lowercase oxide columns listed in *oxides*.
    oxides : list of str
        Oxide column names.  'feo_tot' is treated as FeO.
    basis : float
        Target cation mole sum, typically 100 (default).

    Returns
    -------
    pandas.DataFrame
        Copy of *data* with added cation mole fraction columns and
        column ``'total_cat'``.
    """
    data = data.copy()

    cation_cols = []
    moles = pd.DataFrame(index=data.index)

    for ox in oxides:
        ox_formula = 'FeO' if ox.lower() == 'feo_tot' else ox
        el_sym, n_cat = _parse_cation(ox_formula)
        col = el_sym.capitalize()

        mol_wt = mw('FeO' if ox.lower() == 'feo_tot' else ox_formula)
        moles[col] = data[ox.lower()].clip(lower=0) * n_cat / mol_wt
        cation_cols.append(col)

    total = moles.sum(axis=1)
    total = total.replace(0, np.nan)

    for col in cation_cols:
        data[col] = basis * moles[col] / total
        data.loc[data[col] <= 0, col] = np.nan

    data['total_cat'] = total

    return data


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _count_cations(el_sym, oxide_formula):
    """Number of cation atoms per formula unit (e.g. Al2O3 → 2)."""
    # Match the element symbol at start of formula followed by an optional digit
    pattern = r'^' + re.escape(el_sym) + r'(\d*)'
    m = re.match(pattern, oxide_formula, re.IGNORECASE)
    if m and m.group(1):
        return int(m.group(1))
    return 1


def _parse_cation(oxide_formula):
    """Return (element_symbol, n_cations) from an oxide formula.

    Examples: 'SiO2' → ('Si', 1), 'Al2O3' → ('Al', 2), 'Fe2O3' → ('Fe', 2).
    """
    m = re.match(r'([A-Z][a-z]?)(\d*)', oxide_formula)
    if not m:
        raise ValueError(f"Cannot parse oxide formula: {oxide_formula}")
    sym = m.group(1)
    n   = int(m.group(2)) if m.group(2) else 1
    return sym, n
