"""
Importable pipeline wrapper for the global geochemistry database.

run_pipeline() executes the full load-and-process workflow and returns a
DataFrame.  It is the programmatic equivalent of running
workflows/scripts/gchemload.py.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def run_pipeline(
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

    Parameters
    ----------
    db_version : str
        Database export date tag used to build the default file path.
    normalization : {'anhydrous', 'hydrous', 'none'}
    total_tol : float
        Oxide sum tolerance (100 ± tol wt%).
    ree_scheme : {'reduced', 'full'}
    derived_props : bool
        Run physical property estimates (Vp, Vs, density, Tc, HP).
    hp_method : {'mcdonough', 'rybach'}
    hp_include_rb, hp_include_sm : bool
    hp_est_missing, hp_censor_bdl : bool
    verbose : bool

    Returns
    -------
    pandas.DataFrame
        Fully processed database.
    """
    from src.geochem.fileio.database          import load_database
    from src.geochem.processing.iron          import fefix
    from src.geochem.processing.normalize     import oxide_norm, cation_to_oxide
    from src.geochem.processing.elements      import sumree
    from src.geochem.processing.ages          import age_correction
    from src.geochem.classification.indices   import geochem_index, lambdaree
    from src.geochem.classification.igneous   import tas, graniteclass
    from src.geochem.classification.sedimentary import sedclass
    from src.geochem.classification.metamorphic import metamorphic_class, adjust_origin
    from src.geochem.physprop.seismic         import vpest, vsest
    from src.geochem.physprop.density         import densest
    from src.geochem.physprop.thermal         import tcest
    from src.geochem.physprop.heat            import hpest, hpest_u_th, basalt_liquidus

    def _log(msg):
        if verbose:
            print(msg)

    OXIDES = ['sio2', 'tio2', 'al2o3', 'feo_tot', 'mgo', 'cao', 'na2o', 'k2o', 'p2o5']

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
    _log('Converting Fe to FeO...')
    data = fefix(data)

    # 3. Cation → oxide
    _log('Converting cation data to oxides...')
    data = cation_to_oxide(data)

    # 4. Normalization
    if normalization == 'none':
        _log('Skipping normalization.')
    elif normalization == 'hydrous':
        vol_cols = [c for c in ['loi', 'h2o_plus', 'h2o', 'co2'] if c in data.columns]
        has_vol  = data[vol_cols].notna().any(axis=1)
        data = data[has_vol].reset_index(drop=True)
        _log(f'  {len(data):,} samples with volatile data retained.')
        _log('Normalizing to hydrous conditions...')
        data = oxide_norm(data, oxides=OXIDES, normalization='hydrous', total_tol=total_tol)
    else:
        _log('Normalizing to anhydrous conditions...')
        data = oxide_norm(data, oxides=OXIDES, normalization='anhydrous', total_tol=total_tol)

    # 5. REE sums
    _log(f'Computing REE sums ({ree_scheme} scheme)...')
    data = sumree(data, scheme=ree_scheme)

    # 6. Geochemical indices
    _log('Computing geochemical indices...')
    data = geochem_index(data)
    _log("Fitting REE lambda coefficients (O'Neill 2016)...")
    data = lambdaree(data)

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

    _log('Computing granite classification...')
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

    _log('Inferring protolith origin for metamorphic rocks...')
    data = adjust_origin(data)

    _log('Classifying metamorphic rocks...')
    meta = _metamorphic_mask(data)
    for col in ('met_facies', 'met_texture'):
        data[col] = ''
    if meta.any():
        result = metamorphic_class(data[meta])
        for col in ['met_facies', 'met_texture']:
            if col in result.columns:
                data.loc[result.index, col] = result[col]

    # 8. Physical properties
    if derived_props:
        _log('Computing physical property estimates...')
        data = vpest(data)
        vp_vals = data['p_velocity'].to_numpy(float) if 'p_velocity' in data.columns else None
        if vp_vals is not None:
            data['s_velocity'] = vsest(vp_vals)
        data = densest(data)
        data = tcest(data)
        data = hpest(data, method=hp_method,
                     include_rb=hp_include_rb, include_sm=hp_include_sm,
                     estimate_missing=hp_est_missing, censor_bdl=hp_censor_bdl)
        data = hpest_u_th(data)
        ign = _igneous_mask(data)
        data['basalt_liquidus'] = np.nan
        data.loc[ign, 'basalt_liquidus'] = basalt_liquidus(data[ign])
    else:
        _log('Skipping physical property estimates.')

    _log('Pipeline complete.')
    return data
