"""
Load and process the global whole-rock geochemical database.

Equivalent to the MATLAB gchemload.m workflow.  Designed to be run
interactively section-by-section in the VS Code Python Interactive window
(Shift+Enter on each block), or executed as a whole.

Outputs
-------
data : pd.DataFrame
    Fully processed database with Fe correction, normalization,
    geochemical indices, rock classifications, and physical property
    estimates.

Usage
-----
    Run all at once:
        python scripts/gchemload.py

    Or open in VS Code and run sections with Shift+Enter.
"""

import sys
import pathlib

# Make src/ importable when running from repo root or workflows/scripts/
_root = pathlib.Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(_root))

import numpy as np
import pandas as pd

from src.geochem.fileio.database    import load_database
from src.geochem.processing.iron      import fefix
from src.geochem.processing.normalize import oxide_norm, cation_to_oxide
from src.geochem.processing.elements  import sumree
from src.geochem.processing.ages      import age_correction
from src.geochem.classification.indices    import geochem_index, lambdaree
from src.geochem.classification.igneous    import tas, graniteclass
from src.geochem.classification.sedimentary import sedclass
from src.geochem.classification.metamorphic import metamorphic_class, adjust_origin
from src.geochem.physprop.seismic  import vpest, vsest
from src.geochem.physprop.density  import densest
from src.geochem.physprop.thermal  import tcest
from src.geochem.physprop.heat     import hpest, hpest_u_th, basalt_liquidus


# ============================================================
# Options
# ============================================================

DB_VERSION    = '2024_01_23'
NORMALIZATION = 'anhydrous'   # 'anhydrous', 'hydrous', or 'none'
TOTAL_TOL     = 10            # oxide sum tolerance: 100 ± tol wt%
REE_SCHEME    = 'reduced'     # 'reduced' or 'full'
DERIVED_PROPS = True          # compute physical property estimates

# Heat production options
HP_METHOD       = 'mcdonough'  # 'mcdonough' (default) or 'rybach'
HP_INCLUDE_RB   = False        # include Rb contribution (~0.15% of UCC total)
HP_INCLUDE_SM   = False        # include Sm contribution (~0.04% of UCC total)
HP_EST_MISSING  = True         # predict missing U/Th from Th/U ratio; Rb from K; Sm from Nd
HP_CENSOR_BDL   = True         # apply gaussian censoring to below-detection-limit values

OXIDES = ['sio2', 'tio2', 'al2o3', 'feo_tot', 'mgo', 'cao', 'na2o', 'k2o', 'p2o5']


# ============================================================
# Helpers
# ============================================================

def _group(df, *groups):
    """Boolean mask for rows whose rock_group matches any of *groups* (case-insensitive)."""
    return df['rock_group'].str.lower().isin(g.lower() for g in groups)

def _origin(df, *origins):
    """Boolean mask for rows whose rock_origin matches any of *origins* (case-insensitive)."""
    return df['rock_origin'].str.lower().isin(o.lower() for o in origins)

def _igneous_mask(df):
    """True for igneous rocks and rocks with an igneous protolith."""
    return _group(df, 'igneous') | _origin(df, 'metaigneous', 'metavolcanic', 'metaplutonic')

def _sedimentary_mask(df):
    """True for sedimentary rocks and rocks with a sedimentary protolith."""
    return _group(df, 'sedimentary') | _origin(df, 'metasedimentary')

def _metamorphic_mask(df):
    """True for metamorphic rocks."""
    return _group(df, 'metamorphic')


# ============================================================
# 1. Load database
# ============================================================

print(f'Reading database {DB_VERSION}...')
data = load_database()
print(f'  {len(data):,} samples loaded.\n')

# Drop mineral analyses (not whole-rock)
if 'rock_group' in data.columns:
    data = data[~_group(data, 'mineral')].reset_index(drop=True)
    print(f'  {len(data):,} samples after removing minerals.\n')

# Clean and reconcile age columns
if 'age' in data.columns:
    print('Correcting age data...')
    data = age_correction(data)
    print('  Done.\n')


# ============================================================
# 2. Fe correction
# ============================================================

print('Converting Fe to FeO and computing Fe²⁺/Fe_total...')
data = fefix(data)
print('  Done.\n')


# ============================================================
# 3. Cation → oxide conversion
# ============================================================

print('Converting cation data to oxides where missing...')
data = cation_to_oxide(data)
print('  Done.\n')


# ============================================================
# 4. Oxide normalization
# ============================================================

if NORMALIZATION == 'none':
    print('Skipping normalization.\n')

elif NORMALIZATION == 'hydrous':
    print('Filtering to samples with LOI / volatile data...')
    vol_cols = [c for c in ['loi', 'h2o_plus', 'h2o', 'co2'] if c in data.columns]
    has_vol  = data[vol_cols].notna().any(axis=1)
    data = data[has_vol].reset_index(drop=True)
    print(f'  {len(data):,} samples retained.')
    print(f'Normalizing to hydrous conditions...')
    data = oxide_norm(data, oxides=OXIDES, normalization='hydrous', total_tol=TOTAL_TOL)
    print('  Done.\n')

else:  # anhydrous
    print('Normalizing to anhydrous conditions...')
    data = oxide_norm(data, oxides=OXIDES, normalization='anhydrous', total_tol=TOTAL_TOL)
    print('  Done.\n')


# ============================================================
# 5. REE sums
# ============================================================

print(f'Computing REE sums ({REE_SCHEME} scheme)...')
data = sumree(data, scheme=REE_SCHEME)
print('  Done.\n')


# ============================================================
# 6. Geochemical indices
# ============================================================

print('Computing geochemical indices...')
data = geochem_index(data)
print('  Done.\n')

print('Fitting REE lambda coefficients (O\'Neill 2016)...')
data = lambdaree(data)
n_ree = int(data['nree'].notna().sum())
print(f'  {n_ree:,} samples with REE fits.\n')


# ============================================================
# 7. Rock classification
# ============================================================

# --- TAS (igneous) ---
print('Classifying igneous rocks by TAS...')
ign  = _igneous_mask(data)
rock_category = pd.Series('', index=data.index, dtype=object)
if 'rock_origin' in data.columns:
    rock_category[ign] = data.loc[ign, 'rock_origin'].str.lower()

data['tas_name'] = ''
data.loc[ign, 'tas_name'] = tas(
    data[ign],
    rock_category=rock_category[ign] if ign.any() else None,
)
print(f"  {(data['tas_name'] != '').sum():,} samples classified.\n")

# --- Granite geochemical classification ---
print('Computing granite classification (Frost et al. 2001 + SIA)...')
gc = graniteclass(data[ign])
data['sia_type']     = ''
data['fe_mg_class']  = ''
data['alkali_lime']  = ''
data['alumina_sat']  = ''
data.loc[ign, 'sia_type']    = gc['sia'].values
data.loc[ign, 'fe_mg_class'] = gc['fe_mg'].values
data.loc[ign, 'alkali_lime'] = gc['alkali_lime'].values
data.loc[ign, 'alumina_sat'] = gc['alumina_sat'].values
print('  Done.\n')

# --- Sedimentary classification ---
print('Classifying sedimentary rocks...')
sed = _sedimentary_mask(data)
for col, default in [('sed_name', ''), ('sed_Q', np.nan), ('sed_F', np.nan), ('sed_L', np.nan)]:
    data[col] = default
if sed.any():
    result = sedclass(data[sed])
    for col in ['sed_name', 'sed_Q', 'sed_F', 'sed_L']:
        if col in result.columns:
            data.loc[result.index, col] = result[col]
print(f"  {(data['sed_name'] != '').sum():,} samples classified.\n")

# --- Adjust rock origin for metamorphic samples ---
print('Inferring protolith origin for metamorphic rocks...')
data = adjust_origin(data)
print('  Done.\n')

# --- Metamorphic classification ---
print('Classifying metamorphic rocks (facies + texture)...')
meta = _metamorphic_mask(data)
for col in ('met_facies', 'met_texture'):
    data[col] = ''
if meta.any():
    result = metamorphic_class(data[meta])
    for col in ['met_facies', 'met_texture']:
        if col in result.columns:
            data.loc[result.index, col] = result[col]
print(f"  {(data['met_facies'] != '').sum():,} samples classified.\n")


# ============================================================
# 8. Physical property estimates
# ============================================================

if DERIVED_PROPS:
    print('Computing physical property estimates...')

    print('  P-velocity...')
    data = vpest(data)

    print('  S-velocity...')
    vp_vals = data['p_velocity'].to_numpy(float) if 'p_velocity' in data.columns else None
    if vp_vals is not None:
        data['s_velocity'] = vsest(vp_vals)

    print('  Density...')
    data = densest(data)

    print('  Thermal conductivity...')
    data = tcest(data)

    print(f'  Heat production ({HP_METHOD!r}, est_missing={HP_EST_MISSING}, censor_bdl={HP_CENSOR_BDL})...')
    data = hpest(data, method=HP_METHOD,
                 include_rb=HP_INCLUDE_RB, include_sm=HP_INCLUDE_SM,
                 estimate_missing=HP_EST_MISSING, censor_bdl=HP_CENSOR_BDL)

    print('  Diagnostic U/Th estimates (ratio_th_u, u_est, th_est columns)...')
    data = hpest_u_th(data)

    print('  Basalt liquidus temperature...')
    ign = _igneous_mask(data)
    data['basalt_liquidus'] = np.nan
    data.loc[ign, 'basalt_liquidus'] = basalt_liquidus(data[ign])
    n_liq = int(data['basalt_liquidus'].notna().sum())
    print(f'    {n_liq:,} basaltic samples with liquidus estimate.')

    print('  Done.\n')

else:
    print('Skipping physical property estimates.\n')


# ============================================================
# Summary
# ============================================================

print('=' * 50)
print(f'Database ready: {len(data):,} samples, {len(data.columns)} columns')
if 'rock_group' in data.columns:
    print('\nSamples by rock group:')
    print(data['rock_group'].str.lower().value_counts().to_string())
print('=' * 50)
