"""
Full load-and-process pipeline for the global geochemical database.

``load_geochem()`` is the single entry point for loading the database and
applying the complete processing chain: Fe correction, oxide normalisation,
REE summation, geochemical indices, rock classification, and physical property
estimates.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from global_geochemistry.geochem.fileio.database          import load_database
from global_geochemistry.geochem.processing.iron          import fefix
from global_geochemistry.geochem.processing.normalize     import oxide_norm, cation_to_oxide
from global_geochemistry.geochem.processing.elements      import sumree
from global_geochemistry.geochem.processing.ages          import age_correction
from global_geochemistry.geochem.classification.indices   import geochem_index, lambdaree
from global_geochemistry.geochem.classification.igneous   import tas, graniteclass
from global_geochemistry.geochem.classification.sedimentary import sedclass
from global_geochemistry.geochem.classification.metamorphic import metamorphic_class, adjust_origin
from global_geochemistry.geochem.physprop.seismic         import vpest, vsest
from global_geochemistry.geochem.physprop.density         import densest
from global_geochemistry.geochem.physprop.thermal         import tcest
from global_geochemistry.geochem.physprop.heat            import hpest, hpest_u_th, basalt_liquidus


_OXIDES = ['sio2', 'tio2', 'al2o3', 'feo_tot', 'mgo', 'cao', 'na2o', 'k2o', 'p2o5']


def _group(df, *groups):
    return df['rock_group'].str.lower().isin(g.lower() for g in groups)

def _origin(df, *origins):
    return df['rock_origin'].str.lower().isin(o.lower() for o in origins)

def _igneous_mask(df):
    return _group(df, 'igneous') | _origin(df, 'metaigneous', 'metavolcanic', 'metaplutonic')

def _sedimentary_mask(df):
    return _group(df, 'sedimentary') | _origin(df, 'metasedimentary')

def _metamorphic_mask(df):
    return _group(df, 'metamorphic')


def load_geochem(
    db_version: str = '2024_01_23',
    normalization: str = 'anhydrous',
    total_tol: float = 10,
    ree_scheme: str = 'reduced',
    derived_props: bool = True,
    hp_method: str = 'mcdonough',
    hp_include_rb: bool = False,
    hp_include_sm: bool = False,
    hp_est_missing: bool = True,
    hp_censor_bdl: bool = True,
    verbose: bool = True,
) -> pd.DataFrame:
    """Load and fully process the global geochemical database.

    Applies the complete processing chain in order:

    1. Load raw database CSV
    2. Drop mineral analyses; reconcile age data
    3. Fe correction (fefix)
    4. Cation → oxide conversion
    5. Oxide normalisation (anhydrous / hydrous / none)
    6. REE sums
    7. Geochemical indices + REE λ coefficients
    8. Rock classification (TAS, granite, sedimentary, metamorphic)
    9. Physical property estimates (Vp, Vs, density, Tc, HP) if *derived_props*

    Parameters
    ----------
    db_version : str
        Database export date tag used to build the default CSV path.
    normalization : {'anhydrous', 'hydrous', 'none'}
        Oxide normalisation mode.  ``'hydrous'`` retains only samples with LOI
        or volatile data.
    total_tol : float
        Oxide sum tolerance; samples with total outside ``100 ± total_tol``
        wt% are filtered before normalisation.
    ree_scheme : {'reduced', 'full'}
        Element set used by :func:`~.elements.sumree`.
    derived_props : bool
        Compute physical property estimates (Vp, Vs, density, Tc, HP).
        Set to ``False`` to skip and save time when only geochemistry is needed.
    hp_method : {'mcdonough', 'rybach'}
        Heat production formula.
    hp_include_rb : bool
        Include Rb contribution to heat production (~0.15% of UCC total).
    hp_include_sm : bool
        Include Sm contribution (~0.04% of UCC total).
    hp_est_missing : bool
        Predict missing U/Th from rock-type Th/U ratio; Rb from K; Sm from Nd.
    hp_censor_bdl : bool
        Apply Gaussian censoring to below-detection-limit U, Th, K values.
    verbose : bool
        Print progress messages to stdout.

    Returns
    -------
    pandas.DataFrame
        Fully processed database.

    Examples
    --------
    >>> data = load_geochem()                         # all defaults
    >>> data = load_geochem(derived_props=False)      # skip physprop estimates
    >>> data = load_geochem(normalization='hydrous')  # volatile-bearing samples only
    """
    def _log(msg):
        if verbose:
            print(msg)

    # 1. Load
    _log(f'Reading database {db_version}...')
    data = load_database()
    _log(f'  {len(data):,} samples loaded.')

    if 'rock_group' in data.columns:
        data = data[~_group(data, 'mineral')].reset_index(drop=True)
        _log(f'  {len(data):,} samples after removing minerals.')

    if 'age' in data.columns:
        _log('Correcting age data...')
        data = age_correction(data)

    # 2. Fe correction
    _log('Converting Fe to FeO and computing Fe²⁺/Fe_total...')
    data = fefix(data)

    # 3. Cation → oxide
    _log('Converting cation data to oxides where missing...')
    data = cation_to_oxide(data)

    # 4. Normalisation
    if normalization == 'none':
        _log('Skipping normalisation.')
    elif normalization == 'hydrous':
        vol_cols = [c for c in ['loi', 'h2o_plus', 'h2o', 'co2'] if c in data.columns]
        has_vol  = data[vol_cols].notna().any(axis=1)
        data = data[has_vol].reset_index(drop=True)
        _log(f'  {len(data):,} samples with volatile data retained.')
        _log('Normalising to hydrous conditions...')
        data = oxide_norm(data, oxides=_OXIDES, normalization='hydrous', total_tol=total_tol)
    else:
        _log('Normalising to anhydrous conditions...')
        data = oxide_norm(data, oxides=_OXIDES, normalization='anhydrous', total_tol=total_tol)

    # 5. REE sums
    _log(f'Computing REE sums ({ree_scheme} scheme)...')
    data = sumree(data, scheme=ree_scheme)

    # 6. Geochemical indices
    _log('Computing geochemical indices...')
    data = geochem_index(data)
    _log("Fitting REE lambda coefficients (O'Neill 2016)...")
    data = lambdaree(data)
    if 'nree' in data.columns:
        _log(f"  {int(data['nree'].notna().sum()):,} samples with REE fits.")

    # 7. Rock classification
    _log('Classifying igneous rocks by TAS...')
    ign = _igneous_mask(data)
    rock_category = pd.Series('', index=data.index, dtype=object)
    if 'rock_origin' in data.columns:
        rock_category[ign] = data.loc[ign, 'rock_origin'].str.lower()
    data['tas_name'] = ''
    data.loc[ign, 'tas_name'] = tas(
        data[ign],
        rock_category=rock_category[ign] if ign.any() else None,
    )
    _log(f"  {(data['tas_name'] != '').sum():,} samples classified.")

    _log('Computing granite classification (Frost et al. 2001 + SIA)...')
    gc = graniteclass(data[ign])
    for col in ('sia_type', 'fe_mg_class', 'alkali_lime', 'alumina_sat'):
        data[col] = ''
    data.loc[ign, 'sia_type']    = gc['sia'].values
    data.loc[ign, 'fe_mg_class'] = gc['fe_mg'].values
    data.loc[ign, 'alkali_lime'] = gc['alkali_lime'].values
    data.loc[ign, 'alumina_sat'] = gc['alumina_sat'].values

    _log('Classifying sedimentary rocks...')
    sed = _sedimentary_mask(data)
    for col, default in [('sed_name', ''), ('sed_Q', np.nan), ('sed_F', np.nan), ('sed_L', np.nan)]:
        data[col] = default
    if sed.any():
        result = sedclass(data[sed])
        for col in ['sed_name', 'sed_Q', 'sed_F', 'sed_L']:
            if col in result.columns:
                data.loc[result.index, col] = result[col]
    _log(f"  {(data['sed_name'] != '').sum():,} samples classified.")

    _log('Inferring protolith origin for metamorphic rocks...')
    data = adjust_origin(data)

    _log('Classifying metamorphic rocks (facies + texture)...')
    meta = _metamorphic_mask(data)
    for col in ('met_facies', 'met_texture'):
        data[col] = ''
    if meta.any():
        result = metamorphic_class(data[meta])
        for col in ['met_facies', 'met_texture']:
            if col in result.columns:
                data.loc[result.index, col] = result[col]
    _log(f"  {(data['met_facies'] != '').sum():,} samples classified.")

    # 8. Physical property estimates
    if derived_props:
        _log('Computing physical property estimates...')
        _log('  P-velocity...')
        data = vpest(data)
        _log('  S-velocity...')
        if 'p_velocity' in data.columns:
            data['s_velocity'] = vsest(data['p_velocity'].to_numpy(float))
        _log('  Density...')
        data = densest(data)
        _log('  Thermal conductivity...')
        data = tcest(data)
        _log(f'  Heat production ({hp_method!r}, est_missing={hp_est_missing}, censor_bdl={hp_censor_bdl})...')
        data = hpest(data, method=hp_method,
                     include_rb=hp_include_rb, include_sm=hp_include_sm,
                     estimate_missing=hp_est_missing, censor_bdl=hp_censor_bdl)
        _log('  Diagnostic U/Th estimates...')
        data = hpest_u_th(data)
        _log('  Basalt liquidus temperature...')
        ign = _igneous_mask(data)
        data['basalt_liquidus'] = np.nan
        data.loc[ign, 'basalt_liquidus'] = basalt_liquidus(data[ign])
        if 'basalt_liquidus' in data.columns:
            _log(f"    {int(data['basalt_liquidus'].notna().sum()):,} basaltic samples with liquidus estimate.")
    else:
        _log('Skipping physical property estimates.')

    _log('Done.')
    return data
