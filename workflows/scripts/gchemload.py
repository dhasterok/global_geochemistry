"""
Load and process the global whole-rock geochemical database.

Edit the Options block below, then run:
    python workflows/scripts/gchemload.py

Or open in VS Code and run the whole file with Shift+Enter.
Outputs ``data``, a fully processed pandas DataFrame.
"""

import sys
import pathlib

_root = pathlib.Path(__file__).resolve().parent.parent.parent
if str(_root) not in sys.path:
    sys.path.insert(0, str(_root))

from src.geochem.fileio.pipeline import load_geochem  # noqa: E402


# ============================================================
# Options
# ============================================================

DB_VERSION    = '2024_01_23'
NORMALIZATION = 'anhydrous'   # 'anhydrous', 'hydrous', or 'none'
TOTAL_TOL     = 10            # oxide sum tolerance: 100 ± tol wt%
REE_SCHEME    = 'reduced'     # 'reduced' or 'full'
DERIVED_PROPS = True          # compute physical property estimates

# Heat production options
HP_METHOD      = 'mcdonough'  # 'mcdonough' (default) or 'rybach'
HP_INCLUDE_RB  = False        # include Rb contribution (~0.15% of UCC total)
HP_INCLUDE_SM  = False        # include Sm contribution (~0.04% of UCC total)
HP_EST_MISSING = True         # predict missing U/Th from Th/U ratio; Rb from K; Sm from Nd
HP_CENSOR_BDL  = True         # apply gaussian censoring to below-detection-limit values


# ============================================================
# Run pipeline
# ============================================================

data = load_geochem(
    db_version     = DB_VERSION,
    normalization  = NORMALIZATION,
    total_tol      = TOTAL_TOL,
    ree_scheme     = REE_SCHEME,
    derived_props  = DERIVED_PROPS,
    hp_method      = HP_METHOD,
    hp_include_rb  = HP_INCLUDE_RB,
    hp_include_sm  = HP_INCLUDE_SM,
    hp_est_missing = HP_EST_MISSING,
    hp_censor_bdl  = HP_CENSOR_BDL,
    verbose        = True,
)


# ============================================================
# Summary
# ============================================================

print('=' * 50)
print(f'Database ready: {len(data):,} samples, {len(data.columns)} columns')
if 'rock_group' in data.columns:
    print('\nSamples by rock group:')
    print(data['rock_group'].str.lower().value_counts().to_string())
print('=' * 50)
