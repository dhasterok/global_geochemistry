"""
Chemical classification of sedimentary rocks.

Stage 1: QFL ternary diagram using molar Q (SiO2), F (Al2O3+Fe2O3), L (MgO+CaO).
Stage 2: log(SiO2/Al2O3) vs log(Fe2O3/K2O) for clastic sandstone sub-types.

References:
    Herron, M.M. (1988) Geochemical classification of terrigenous sands and
    shales from core or log data.  Journal of Sedimentary Research 58(5),
    820–829.  (log-space sandstone discriminant)

    Folk, R.L. (1980) Petrology of Sedimentary Rocks.  Hemphill Publishing,
    Austin.  (QFL framework, slightly modified here for geochemistry).
"""

import numpy as np
import pandas as pd
from matplotlib.path import Path

from src.utils.molecular import MolecularWeightCalculator

_mwc = MolecularWeightCalculator()

def _mw(f): return _mwc.molecular_weight(f)

_MW = {ox: _mw(ox) for ox in ['SiO2', 'Al2O3', 'Fe2O3', 'MgO', 'CaO', 'K2O', 'FeO']}

# Conversion factor: FeO_tot → Fe2O3_tot (wt ratio)
_FE_CF = _MW['Fe2O3'] / (2 * _MW['FeO'])

# ---------------------------------------------------------------------------
# QFL polygon definitions (Q, F, L coordinates, each 0–100%)
# Source: load_sedgons.m
# ---------------------------------------------------------------------------
_QFL_POLYGONS = [
    (np.array([[0,0],[0,10],[30,35],[40,30],[58,12],[55,0]]),
     'carbonate'),
    (np.array([[100,0],[90,10],[90,0]]),
     'quartzite'),
    (np.array([[90,10],[90,0],[55,0],[58,12],[83,17]]),
     'psammite'),
    (np.array([[83,17],[58,12],[40,30],[30,35],[30,60],[30,70]]),
     'pelite'),
    (np.array([[0,100],[0,95],[30,60],[30,70]]),
     'laterite'),
    (np.array([[0,10],[30,35],[30,60],[0,95]]),
     'oxide'),
]
_QFL_PATHS = [Path(v) for v, _ in _QFL_POLYGONS]

# ---------------------------------------------------------------------------
# Log-space sandstone polygon definitions (log10 SiO2/Al2O3, log10 Fe2O3/K2O)
# Source: load_sedgons2.m
# ---------------------------------------------------------------------------
_SAND_POLYGONS = [
    (np.array([[-2,0.6],[0.71,0.6],[0.71,5],[-2,5]]),
     'iron-rich shale'),
    (np.array([[0.71,0.6],[1.728,0.6],[2.08,5],[0.71,5]]),
     'iron-rich sand'),
    (np.array([[1.52,-2],[2.08,5],[4,5],[4,-2]]),
     'quartz arenite'),
    (np.array([[1.684,0.05],[1.092,0.05],[1.1375,0.6],[1.728,0.6]]),
     'sublitharenite'),
    (np.array([[0.9124,-2],[1.092,0.05],[1.684,0.05],[1.52,-2]]),
     'subarkose'),
    (np.array([[0.3061,-2],[0.7665,0.05],[1.092,0.05],[0.9124,-2]]),
     'arkose'),
    (np.array([[0.7665,0.05],[0.89,0.6],[1.1375,0.6],[1.092,0.05]]),
     'litharenite'),
    (np.array([[0.1157,-2],[0.71,0.6],[0.89,0.6],[0.3061,-2]]),
     'wacke'),
    (np.array([[0.1157,-2],[0.71,0.6],[-2,0.6],[-2,-2]]),
     'shale'),
]
_SAND_PATHS = [Path(v) for v, _ in _SAND_POLYGONS]


def sedclass(data):
    """Classify sedimentary rocks from major-oxide composition.

    Two-stage classification:
    1. QFL ternary (SiO2, Al2O3+Fe2O3, MgO+CaO) → broad lithology.
    2. log(SiO2/Al2O3) vs log(Fe2O3/K2O) → sandstone sub-type for clastic
       samples (pelite / psammite) not already assigned a carbonate name.

    Parameters
    ----------
    data : pandas.DataFrame
        Must contain 'sio2', 'al2o3', 'feo_tot', 'mgo', 'cao', 'k2o'.

    Returns
    -------
    pandas.DataFrame
        Input data with three new columns:
        - 'sed_name': rock name string
        - 'sed_Q', 'sed_F', 'sed_L': molar QFL coordinates (0–100%)
    """
    n = len(data)

    def col(name):
        arr = data[name].to_numpy(float) if name in data.columns else np.zeros(n)
        return np.where(arr < 0, np.nan, arr)

    sio2   = col('sio2')
    al2o3  = col('al2o3')
    feo_t  = col('feo_tot')
    mgo    = col('mgo')
    cao    = col('cao')
    k2o    = col('k2o')

    fe2o3 = feo_t * _FE_CF

    # QFL mole fractions
    Q = sio2  / _MW['SiO2']
    F = al2o3 / _MW['Al2O3'] + fe2o3 / _MW['Fe2O3']
    L = mgo   / _MW['MgO']   + cao   / _MW['CaO']

    QFL_total = Q + F + L
    with np.errstate(invalid='ignore', divide='ignore'):
        Qn = np.where(QFL_total > 0, Q / QFL_total * 100, np.nan)
        Fn = np.where(QFL_total > 0, F / QFL_total * 100, np.nan)
        Ln = np.where(QFL_total > 0, L / QFL_total * 100, np.nan)

    out = data.copy()
    out['sed_Q'] = Qn
    out['sed_F'] = Fn
    out['sed_L'] = Ln

    names = np.full(n, '', dtype=object)

    # Stage 1: QFL ternary classification
    valid_qfl = ~(np.isnan(Qn) | np.isnan(Fn))
    pts_qfl = np.column_stack([Qn[valid_qfl], Fn[valid_qfl]])

    for path, (_, label) in zip(_QFL_PATHS, _QFL_POLYGONS):
        inside = np.zeros(n, dtype=bool)
        if pts_qfl.shape[0]:
            inside[valid_qfl] = path.contains_points(pts_qfl)

        if label == 'carbonate':
            # Dolomite vs limestone based on Mg/(Mg+Ca) mole ratio
            n_mg = mgo / _MW['MgO']
            n_ca = cao / _MW['CaO']
            with np.errstate(invalid='ignore', divide='ignore'):
                mg_frac = np.where((n_mg + n_ca) > 0, n_mg / (n_mg + n_ca), 0.0)
            names[inside & (mg_frac > 0.4)] = 'dolomite'
            names[inside & (mg_frac <= 0.4)] = 'limestone'

        elif label == 'laterite':
            # Bauxite (Al dominant) vs laterite (Fe dominant)
            n_al = al2o3 / _MW['Al2O3']
            n_fe = fe2o3 / _MW['Fe2O3']
            names[inside & (n_al > n_fe)] = 'bauxite'
            names[inside & (n_al <= n_fe)] = 'laterite'

        else:
            names[inside] = label

    # Stage 2: log-space sandstone sub-classification
    # Applied to pelite/psammite samples with required positive values
    clastic = np.isin(names, ['pelite', 'psammite'])
    with np.errstate(invalid='ignore', divide='ignore'):
        log_sa = np.log10(sio2 / al2o3)
        log_fk = np.log10(fe2o3 / k2o)

    valid_sand = clastic & np.isfinite(log_sa) & np.isfinite(log_fk)
    if valid_sand.any():
        pts_sand = np.column_stack([log_sa[valid_sand], log_fk[valid_sand]])
        for path, (_, label) in zip(_SAND_PATHS, _SAND_POLYGONS):
            inside = np.zeros(n, dtype=bool)
            inside[valid_sand] = path.contains_points(pts_sand)
            names[inside] = label

    out['sed_name'] = names
    return out
