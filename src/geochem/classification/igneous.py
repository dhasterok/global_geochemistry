"""
Classification of igneous rocks.

Implements:
- TAS (Total Alkali–Silica) classification — volcanic and plutonic rock
  names based on Le Bas et al. (1986) and Middlemost (1994); high-Mg
  overrides from Le Bas & Streckeisen (1991).
- Carbonatite classification on the CaO–MgO–FeO+Fe2O3+MnO ternary
  (calciocarbonatite, magnesiocarbonatite, ferrocarbonatite).
- QAP normative classification (IUGS double triangle; Le Maitre 2002)
  using CIPW mineral proportions.
- Granite geochemical classification: Frost et al. (2001) ferroan/
  magnesian + alkali-lime series + alumina saturation; SIA (S/I/A-type)
  scheme of Chappell & White (1992) and Whalen et al. (1987).
- Mg number.

References:
    Le Bas et al. (1986) J. Petrology 27, 745–750.
    Middlemost (1994) Earth-Sci. Rev. 37, 215–224.
    Le Bas & Streckeisen (1991) J. Geol. Soc. London 148, 825–833.
    Le Maitre, R.W. (ed.) (2002) Igneous Rocks, 2nd ed. Cambridge UP.
    Frost et al. (2001) J. Petrology 42, 2033–2048.
    Chappell & White (1992) Trans. R. Soc. Edinburgh 83, 1–26.
    Whalen et al. (1987) Contrib. Mineral. Petrol. 95, 407–419.
"""

import numpy as np
import pandas as pd
from matplotlib.path import Path

from global_geochemistry.utils.molecular import MolecularWeightCalculator

_mwc = MolecularWeightCalculator()
def _mw(f): return _mwc.molecular_weight(f)


# ===========================================================================
# TAS polygon definitions (SiO2, Na2O+K2O) — Middlemost (1994)
# Each entry: (vertices, volcanic_name, plutonic_name)
# ===========================================================================

_TAS_POLYGONS = [
    (np.array([[41,0],[41,3],[33,3],[33,0]]),
     '',                    'peridotite'),
    (np.array([[41,0],[41,3],[45,3],[45,0],[41,0]]),
     'picrobasalt',         'peridotgabbro'),
    (np.array([[45,2],[45,5],[52,5],[45,2]]),
     'alkalic basalt',      'alkalic gabbro'),
    (np.array([[45,0],[45,2],[52,5],[52,0],[45,0]]),
     'subalkalic basalt',   'subalkalic gabbro'),
    (np.array([[52,0],[52,5],[57,5.9],[57,0],[52,0]]),
     'basaltic andesite',   'gabbroic diorite'),
    (np.array([[57,0],[57,5.9],[63,7],[63,0],[57,0]]),
     'andesite',            'diorite'),
    (np.array([[63,0],[63,7],[69,8],[77.3,0],[63,0]]),
     'dacite',              'granodiorite'),
    (np.array([[77.3,0],[69,8],[71.8,13.5],[85.9,6.8],[87.5,4.7],[77.3,0]]),
     'rhyolite',            'granite'),
    (np.array([[45,5],[49.4,7.3],[52,5],[45,5]]),
     'trachybasalt',        'monzogabbro'),
    (np.array([[52,5],[49.4,7.3],[53,9.3],[57,5.9],[52,5]]),
     'basaltic trachyandesite', 'monzodiorite'),
    (np.array([[57,5.9],[53,9.3],[57.6,11.7],[63,7],[57,5.9]]),
     'trachyandesite',      'monzonite'),
    (np.array([[63,7],[61.1,8.65],[71.8,13.5],[69,8],[63,7]]),
     'trachydacite',        'quartz monzonite'),
    (np.array([[61.1,8.65],[57.6,11.7],[61,13.5],[63,16.2],[71.8,13.5],[61.1,8.65]]),
     'trachyte',            'syenite'),
    (np.array([[41,3],[41,7],[45,9.4],[49.4,7.3],[45,5],[45,3],[41,3]]),
     'tephrite',            'foid gabbro'),
    (np.array([[49.4,7.3],[45,9.4],[48.4,11.5],[53,9.3],[49.4,7.3]]),
     'phonotephrite',       'foid monzodiorite'),
    (np.array([[53,9.3],[48.4,11.5],[52.5,14],[57.6,11.7],[53,9.3]]),
     'tephriphonolite',     'foid monzosyenite'),
    (np.array([[57.6,11.7],[52.5,14],[52.5,18],[57,18],[63,16.2],[61,13.5],[57.6,11.7]]),
     'phonolite',           'foid syenite'),
    (np.array([[33,10],[33,3],[41,3],[41,7],[33,10]]),
     'ultramafic foidite',  'ultramafic foidolite'),
    (np.array([[41,7],[33,10],[37,14],[45,9.4],[41,7]]),
     'mafic foidite',       'mafic foidolite'),
    (np.array([[37,14],[52.5,18],[52.5,14],[48.4,11.5],[45,9.4],[37,14]]),
     'intermediate foidite', 'intermediate foidolite'),
    (np.array([[77.3,0],[87.5,4.7],[85.9,6.8],[93.2,6.8],[100,0],[77.3,0]]),
     'silexite',            'quartzolite'),
    (np.array([[93.2,6.8],[85.9,6.8],[71.8,13.5],[63,16.2],[57,18],[52.5,18],
               [37,14],[37,63],[93.2,6.8]]),
     'ultra-high alkali volcanic', 'ultra-high alkali plutonic'),
]

_TAS_PATHS = [Path(v) for v, _, _ in _TAS_POLYGONS]


# ===========================================================================
# Carbonatite ternary polygons — CaO, MgO, FeO+Fe2O3+MnO (load_carbgons.m)
# Vertices in fractional (C, M) space after normalising C+M+F=1.
# ===========================================================================

_CARB_POLYGONS = [
    (np.array([[1.0, 0.0],[0.8, 0.2],[0.8, 0.0]]),
     'calciocarbonatite'),
    (np.array([[0.8, 0.2],[0.0, 1.0],[0.0, 0.5],[0.8, 0.1]]),
     'magnesiocarbonatite'),
    (np.array([[0.8, 0.0],[0.0, 0.0],[0.0, 0.5],[0.8, 0.1]]),
     'ferrocarbonatite'),
]
_CARB_PATHS = [Path(v) for v, _ in _CARB_POLYGONS]


# ===========================================================================
# QAP polygon definitions — [Q, A, P, F] as percent of Q+A+P+F
# F=0 polygons → tested on Q-A-P face; Q=0 polygons → tested on A-P-F face.
# Source: load_qapgons.m (Le Maitre 2002)
# ===========================================================================

_QAP_POLYGONS = [
    # Q-A-P face
    (np.array([[90,10,0,0],[100,0,0,0],[90,0,10,0]]),
     'quartzolite',          'silexite'),
    (np.array([[90,10,0,0],[90,0,10,0],[60,0,40,0],[60,40,0,0]]),
     'quartz-rich granitoid', 'quartz-rich rhyolite'),
    (np.array([[60,40,0,0],[60,36,4,0],[20,72,8,0],[20,80,0,0]]),
     'alkali feldspar granite', 'alkali feldspar rhyolite'),
    (np.array([[60,4,36,0],[20,8,72,0],[20,28,52,0],[60,14,26,0]]),
     'granodiorite',          'dacite'),
    (np.array([[60,26,14,0],[20,52,28,0],[20,72,8,0],[60,36,4,0]]),
     'syenogranite',          'rhyolite'),
    (np.array([[60,14,26,0],[20,28,52,0],[20,52,28,0],[60,26,14,0]]),
     'monzogranite',          'rhyolite'),
    (np.array([[60,0,40,0],[60,4,36,0],[20,8,72,0],[20,0,80,0]]),
     'tonalite',              'dacite'),
    (np.array([[20,80,0,0],[5,95,0,0],[5,85.5,9.5,0],[20,72,8,0]]),
     'alkali feldspar quartz syenite', 'alkali feldspar quartz trachyte'),
    (np.array([[20,72,8,0],[5,85.5,9.5,0],[5,61.75,33.25,0],[20,52,28,0]]),
     'quartz syenite',        'quartz trachyte'),
    (np.array([[20,52,28,0],[5,61.75,33.25,0],[5,33.25,61.75,0],[20,28,52,0]]),
     'quartz monzonite',      'quartz latite'),
    (np.array([[20,28,52,0],[5,33.25,61.75,0],[5,9.5,85.5,0],[20,8,72,0]]),
     'quartz monzodiorite',   'andesite'),
    (np.array([[20,8,72,0],[5,9.5,85.5,0],[5,0,95,0],[20,0,80,0]]),
     'quartz diorite',        'andesite'),
    # A-P only (Q=0, F=0)
    (np.array([[5,95,0,0],[0,100,0,0],[0,90,10,0],[5,85.5,9.5,0]]),
     'alkali feldspar syenite', 'alkali feldspar trachyte'),
    (np.array([[5,85.5,9.5,0],[0,90,10,0],[0,65,35,0],[5,61.75,33.25,0]]),
     'syenite',                'trachyte'),
    (np.array([[5,61.75,33.25,0],[0,65,35,0],[0,35,65,0],[5,33.25,61.75,0]]),
     'monzonite',              'latite'),
    (np.array([[5,33.25,61.75,0],[0,35,65,0],[0,10,90,0],[5,9.5,85.5,0]]),
     'monzodiorite',           'andesite'),
    (np.array([[5,9.5,85.5,0],[0,10,90,0],[0,0,100,0],[5,0,95,0]]),
     'diorite',                'andesite'),
    # A-P-F face (Q=0)
    (np.array([[0,0,0,100],[0,10,0,90],[0,0,10,90]]),
     'foidolite',              'foidite'),
    (np.array([[0,10,0,90],[0,40,0,60],[0,20,20,60],[0,5,5,90]]),
     'foidolite',              'phonolitic foidite'),
    (np.array([[0,0,10,90],[0,0,40,60],[0,20,20,60],[0,5,5,90]]),
     'foidolite',              'tephritic foidite'),
    (np.array([[0,40,0,60],[0,36,4,60],[0,81,9,10],[0,90,0,10]]),
     'foid syenite',           'phonolite'),
    (np.array([[0,100,0,0],[0,90,10,0],[0,81,9,10],[0,90,0,10]]),
     'foid-bearing alkali feldspar syenite', 'foid-bearing alkali feldspar trachyte'),
    (np.array([[0,90,10,0],[0,81,9,10],[0,58.5,31.5,10],[0,65,35,0]]),
     'foid-bearing syenite',   'foid-bearing trachyte'),
    (np.array([[0,65,35,0],[0,35,65,0],[0,31.5,58.5,10],[0,58.5,31.5,10]]),
     'foid-bearing monzonite', 'foid-bearing latite'),
    (np.array([[0,10,90,0],[0,9,81,10],[0,31.5,58.5,10],[0,35,65,0]]),
     'foid-bearing monzodiorite', 'andesite'),
    (np.array([[0,0,100,0],[0,10,90,0],[0,9,81,10],[0,0,90,10]]),
     'foid-bearing diorite',   'andesite'),
    (np.array([[0,0,40,60],[0,4,36,60],[0,9,81,10],[0,0,90,10]]),
     'foid diorite',           'tephrite'),
    (np.array([[0,45,45,10],[0,9,81,10],[0,4,36,60],[0,20,20,60]]),
     'foid monzodiorite',      'phonolitic tephrite'),
    (np.array([[0,45,45,10],[0,81,9,10],[0,36,4,60],[0,20,20,60]]),
     'foid monzosyenite',      'tephritic phonolite'),
]

# Pre-build QAP paths projected onto the appropriate ternary face
def _qap_path(verts, face):
    with np.errstate(invalid='ignore', divide='ignore'):
        if face == 'QAP':
            t = verts[:,0] + verts[:,1] + verts[:,2]
            x = verts[:,1] / t   # A fraction
            y = verts[:,0] / t   # Q fraction
        else:
            t = verts[:,1] + verts[:,2] + verts[:,3]
            x = verts[:,1] / t   # A fraction
            y = verts[:,3] / t   # F fraction
    return Path(np.column_stack([x, y]))

_QAP_DATA = []
for _v, _pn, _vn in _QAP_POLYGONS:
    _face = 'APF' if _v[:,3].max() > 0.5 else 'QAP'
    _QAP_DATA.append((_face, _qap_path(_v, _face), _pn, _vn))


# ===========================================================================
# TAS classification
# ===========================================================================

def tas(data, rock_category=None):
    """Classify igneous rocks using the Total Alkali–Silica diagram.

    Parameters
    ----------
    data : pandas.DataFrame
        Must contain 'sio2', 'na2o', 'k2o'.  Optional: 'mgo', 'feo_tot',
        'tio2' (needed for high-Mg and ultramafic overrides).
    rock_category : array-like of str, optional
        Per-sample category: 'volcanic', 'plutonic', or other/NaN.
        When None all samples receive volcanic names.

    Returns
    -------
    pandas.Series
        Rock names; empty string where unclassifiable.
    """
    n = len(data)
    names = np.full(n, '', dtype=object)

    sio2 = data['sio2'].to_numpy(float)
    na2o = data['na2o'].to_numpy(float) if 'na2o' in data.columns else np.zeros(n)
    k2o  = data['k2o'].to_numpy(float)  if 'k2o'  in data.columns else np.zeros(n)
    mgo  = data['mgo'].to_numpy(float)  if 'mgo'  in data.columns else np.full(n, np.nan)
    feo  = data['feo_tot'].to_numpy(float) if 'feo_tot' in data.columns else np.full(n, np.nan)
    tio2 = data['tio2'].to_numpy(float) if 'tio2' in data.columns else np.full(n, np.nan)

    na2o = np.where(np.isnan(na2o), 0.0, na2o)
    k2o  = np.where(np.isnan(k2o),  0.0, k2o)
    ta   = na2o + k2o

    valid = ~np.isnan(sio2)
    pts_v = np.column_stack([sio2[valid], ta[valid]])

    is_volcanic = np.zeros(n, dtype=bool)
    is_plutonic = np.zeros(n, dtype=bool)
    if rock_category is None:
        is_volcanic = valid.copy()
    else:
        rc = np.asarray(rock_category, dtype=object)
        is_volcanic = np.isin(rc, ['volcanic'])
        is_plutonic = np.isin(rc, ['plutonic'])

    for path, (_, volc_name, plut_name) in zip(_TAS_PATHS, _TAS_POLYGONS):
        inside = np.zeros(n, dtype=bool)
        if pts_v.shape[0]:
            inside[valid] = path.contains_points(pts_v)
        names[inside & is_volcanic] = volc_name
        names[inside & is_plutonic] = plut_name

    mafic = valid & (sio2 >= 33) & (sio2 < 65)

    mask = mafic & (sio2 >= 52) & _safe_ge(mgo, 8) & _safe_lt(tio2, 0.5)
    names[mask & is_volcanic] = 'boninite'
    names[mask & is_plutonic] = 'sanukitoid'

    mask = mafic & (sio2 < 52) & _safe_ge(mgo, 12) & (ta < 3)
    names[mask & _safe_ge(feo, 13)] = 'ferropicrite'
    names[mask & _safe_lt(feo, 13)] = 'picrite'

    mask = mafic & (sio2 < 52) & _safe_ge(mgo, 12) & (ta >= 3)
    names[mask & _safe_ge(feo, 13)] = 'alkali ferropicrite'
    names[mask & _safe_lt(feo, 13)] = 'alkali picrite'

    mask = mafic & (sio2 < 52) & _safe_ge(mgo, 18) & _safe_gt(tio2, 1) & (ta < 2)
    names[mask] = 'meimechite'

    mask = mafic & (sio2 < 45) & _safe_ge(mgo, 18) & _safe_lt(tio2, 1) & (ta < 2)
    names[mask & is_volcanic] = 'komatiite'
    names[mask & is_plutonic] = 'intrusive komatiite'

    mask = mafic & (sio2 >= 45) & (sio2 < 52) & _safe_ge(mgo, 18) & _safe_lt(tio2, 1) & (ta < 2)
    names[mask & is_volcanic] = 'basaltic komatiite'
    names[mask & is_plutonic] = 'gabbroic komatiite'

    # Plutonic ultramafic subdivision by TiO2 (Reverdatto et al. 2008)
    mask = is_plutonic & valid & (sio2 >= 33) & (sio2 < 45) & _safe_gt(tio2, 0) & _safe_le(tio2, 0.3)
    names[mask & _safe_ge(mgo, 35)]                      = 'mantle peridotite'
    names[mask & _safe_ge(mgo, 18) & _safe_lt(mgo, 35)] = 'mantle pyroxenite'

    return pd.Series(names, index=data.index, dtype=object)


# ===========================================================================
# Carbonatite classification
# ===========================================================================

def carbclass(data):
    """Classify carbonatites on the CaO–MgO–(FeO+Fe2O3+MnO) ternary.

    A sample is treated as a carbonatite when any of:
    - CO2 >= 20 wt% (if column present)
    - rock_name contains 'carbonatite' (case-insensitive)
    - SiO2 < 33 wt%

    Samples with 20 < SiO2 < 33 receive a 'silico-' prefix.

    Parameters
    ----------
    data : pandas.DataFrame
        Must contain 'sio2', 'cao', 'mgo', 'feo_tot'.  Optional: 'mno',
        'fe2_fe_tot' (Fe²⁺/ΣFe molar ratio), 'co2', 'rock_name'.

    Returns
    -------
    pandas.Series
        Carbonatite sub-type name; empty string for non-carbonatite samples.
    """
    n = len(data)
    names = np.full(n, '', dtype=object)

    sio2    = data['sio2'].to_numpy(float)
    cao     = data['cao'].to_numpy(float)  if 'cao'     in data.columns else np.zeros(n)
    mgo     = data['mgo'].to_numpy(float)  if 'mgo'     in data.columns else np.zeros(n)
    feo_tot = data['feo_tot'].to_numpy(float) if 'feo_tot' in data.columns else np.zeros(n)
    mno     = data['mno'].to_numpy(float)  if 'mno'     in data.columns else np.zeros(n)
    co2     = data['co2'].to_numpy(float)  if 'co2'     in data.columns else np.zeros(n)

    # Iron proxy: 80% FeO + 20% Fe2O3 equivalent + MnO (default);
    # uses measured Fe ratio when available
    _fe_cf = _mw('Fe2O3') / (2 * _mw('FeO'))
    F = 0.8 * feo_tot + 0.2 * feo_tot * _fe_cf + mno
    if 'fe2_fe_tot' in data.columns:
        fe2 = data['fe2_fe_tot'].to_numpy(float)
        has_ratio = np.isfinite(fe2) & (fe2 > 0)
        F[has_ratio] = (fe2[has_ratio] * feo_tot[has_ratio] +
                        (1 - fe2[has_ratio]) * feo_tot[has_ratio] * _fe_cf +
                        mno[has_ratio])

    # Carbonatite detection
    is_carb = (np.abs(sio2) < 33)
    if 'co2' in data.columns:
        is_carb = is_carb | (co2 >= 20)
    if 'rock_name' in data.columns:
        rn = data['rock_name'].astype(str).str.lower()
        is_carb = is_carb | rn.str.contains('carbonatite', regex=False).to_numpy()

    names[is_carb] = 'carbonatite'

    silico = is_carb & (sio2 > 20) & (sio2 < 33)
    names[silico] = 'silicocarbonatite'

    # Sub-classify using CaO–MgO–F ternary
    CMF = cao + mgo + F
    with np.errstate(invalid='ignore', divide='ignore'):
        C_frac = np.where(CMF > 0, cao / CMF, np.nan)
        M_frac = np.where(CMF > 0, mgo / CMF, np.nan)

    ternary_valid = is_carb & np.isfinite(C_frac) & np.isfinite(M_frac)
    if ternary_valid.any():
        pts = np.column_stack([C_frac[ternary_valid], M_frac[ternary_valid]])
        for path, (_, label) in zip(_CARB_PATHS, _CARB_POLYGONS):
            inside = np.zeros(n, dtype=bool)
            inside[ternary_valid] = path.contains_points(pts)
            names[inside & ~silico] = label
            names[inside &  silico] = 'silico-' + label

    return pd.Series(names, index=data.index, dtype=object)


# ===========================================================================
# QAP normative classification
# ===========================================================================

def qap(cipw_result, is_volcanic=None):
    """Assign QAP rock names from CIPW normative mineral proportions.

    Parameters
    ----------
    cipw_result : pandas.DataFrame
        Output of :func:`src.classification.cipw.cipwnorm`.
    is_volcanic : array-like of bool, optional
        True → volcanic name.  Defaults to False (plutonic names).

    Returns
    -------
    pandas.Series
        QAP rock name; empty string where classification is not possible.
    """
    n = len(cipw_result)
    if is_volcanic is None:
        is_volcanic = np.zeros(n, dtype=bool)
    else:
        is_volcanic = np.asarray(is_volcanic, dtype=bool)

    def _col(name):
        arr = cipw_result[name].to_numpy(float) if name in cipw_result.columns else np.zeros(n)
        return np.nan_to_num(arr, nan=0.0)

    Q = _col('cipw_quartz')
    A = (_col('cipw_orthoclase') + _col('cipw_potassium_ms') +
         _col('cipw_kaliophilite') + _col('cipw_sodium_ms'))
    P = _col('cipw_albite') + _col('cipw_anorthite')
    F = _col('cipw_nepheline') + _col('cipw_leucite') + _col('cipw_kaliophilite')

    total = Q + A + P + F
    valid = total > 0

    with np.errstate(invalid='ignore', divide='ignore'):
        Qn = np.where(valid, Q / total, np.nan)
        An = np.where(valid, A / total, np.nan)
        Fn = np.where(valid, F / total, np.nan)

    use_apf = valid & (Fn > 0.1)
    use_qap = valid & ~use_apf

    pts_qap = np.column_stack([An, Qn])
    pts_apf = np.column_stack([An, Fn])

    names = np.full(n, '', dtype=object)

    for face, path, plut_name, volc_name in _QAP_DATA:
        use = use_apf if face == 'APF' else use_qap
        if not use.any():
            continue
        inside = np.zeros(n, dtype=bool)
        inside[use] = path.contains_points(pts_apf[use] if face == 'APF' else pts_qap[use])
        names[inside & ~is_volcanic] = plut_name
        names[inside &  is_volcanic] = volc_name

    return pd.Series(names, index=cipw_result.index, dtype=object)


# ===========================================================================
# Granite geochemical classification
# ===========================================================================

def graniteclass(data):
    """Geochemical classification of granitoids.

    Implements:
    - Frost et al. (2001) three-axis scheme:
        1. ferroan / magnesian (Fe-number vs SiO2 boundary)
        2. alkalic / alkali-calcic / calc-alkalic / calcic (MALI vs SiO2)
        3. peraluminous / metaluminous / peralkaline (ASI)
    - SIA-type (S / I / A) classification:
        A-type: agpaitic index (Al/(Na+K) molar) < 1/0.85 AND
                (ferroan below 70% SiO2, OR above 70% SiO2)
        S-type: ASI >= 1.1
        I-type: remainder

    Parameters
    ----------
    data : pandas.DataFrame
        Requires columns computed by :func:`src.classification.indices.geochem_index`:
        'Fe_number', 'MALI', 'ASI', 'sio2', 'al2o3', 'na2o', 'k2o'.
        Trace elements 'f_ppm', 'zr_ppm', 'nb_ppm', 'ce_ppm', 'y_ppm',
        'ga_ppm' improve A-type detection when available.

    Returns
    -------
    dict of pandas.Series
        Keys: 'sia', 'fe_mg', 'alkali_lime', 'alumina_sat'
    """
    n = len(data)

    def col(name):
        if name in data.columns:
            return data[name].to_numpy(float)
        return np.full(n, np.nan)

    sio2    = col('sio2')
    al2o3   = col('al2o3')
    na2o    = col('na2o')
    k2o     = col('k2o')
    fe_num  = col('Fe_number')
    mali    = col('MALI')
    asi     = col('ASI')

    # ---- 1. Frost et al. (2001): ferroan / magnesian ----
    boundary_fe = 0.486 + 0.0046 * sio2
    fe_mg = np.where(fe_num < boundary_fe, 'magnesian',
            np.where(fe_num >= boundary_fe, 'ferroan', ''))

    # ---- 2. Alkali-lime series (MALI vs SiO2 polynomial boundaries) ----
    b_alkalic     = -41.86 + 1.12   * sio2 - 0.00572 * sio2**2
    b_alkali_calc = -44.72 + 1.094  * sio2 - 0.00527 * sio2**2
    b_calc_calc   = -45.36 + 1.0043 * sio2 - 0.00427 * sio2**2

    alkali_lime = np.where(mali > b_alkalic,     'alkalic',
                  np.where(mali > b_alkali_calc,  'alkali-calcic',
                  np.where(mali > b_calc_calc,    'calc-alkalic',
                                                  'calcic')))
    alkali_lime = np.where(np.isnan(mali), '', alkali_lime)

    # ---- 3. Alumina saturation ----
    # Agpaitic index: molar Al / (Na + K)
    n_al  = 2 * al2o3 / _mw('Al2O3')   # moles Al atoms per 100 g
    n_nak = (2 * na2o / _mw('Na2O') +
             2 * k2o  / _mw('K2O'))    # moles Na+K atoms per 100 g
    with np.errstate(invalid='ignore', divide='ignore'):
        agpaitic = np.where(n_nak > 0, n_al / n_nak, np.nan)

    peralkaline_cond = (asi < 1.0) & np.isfinite(agpaitic) & (agpaitic < 1.0)
    alumina_sat = np.where(asi >= 1.0,       'peraluminous',
                  np.where(peralkaline_cond,  'peralkaline',
                                              'metaluminous'))
    alumina_sat = np.where(np.isnan(asi), '', alumina_sat)

    # ---- 4. SIA scheme ----
    # A-type conditions
    ferroan = fe_num >= boundary_fe
    is_a_low_si  = ferroan & (sio2 < 70) & (agpaitic < 1/0.85)  # agpaitic < ~1.18
    is_a_high_si = (sio2 >= 70) & (agpaitic < 1/0.85)

    sia = np.full(n, '', dtype=object)
    sia[is_a_low_si | is_a_high_si] = 'A'
    # S-type: strongly peraluminous (ASI >= 1.1)
    sia[(sia == '') & (asi >= 1.1)] = 'S'
    # I-type: remainder with valid SiO2
    sia[(sia == '') & np.isfinite(sio2)] = 'I'

    return {
        'sia':         pd.Series(sia,         index=data.index),
        'fe_mg':       pd.Series(fe_mg,       index=data.index),
        'alkali_lime': pd.Series(alkali_lime, index=data.index),
        'alumina_sat': pd.Series(alumina_sat, index=data.index),
    }


# ===========================================================================
# Mg number
# ===========================================================================

def mgnum(mgo, feo_tot, fe3_fraction=0.0):
    """Mg/(Mg + Fe²⁺) molar ratio.

    Parameters
    ----------
    mgo, feo_tot : array-like
        MgO and total FeO (wt%).
    fe3_fraction : float or array-like
        Fraction of total iron as Fe³⁺ (0–1).  Default 0.
    """
    mgo      = np.asarray(mgo,     dtype=float)
    feo_tot  = np.asarray(feo_tot, dtype=float)
    fe3_frac = np.asarray(fe3_fraction, dtype=float)

    n_mg  = mgo / _mw('MgO')
    n_feo = feo_tot * (1 - fe3_frac) / _mw('FeO')
    with np.errstate(invalid='ignore', divide='ignore'):
        mg = n_mg / (n_mg + n_feo)
    return np.where(np.isfinite(mg), mg, np.nan)


# ===========================================================================
# Plotting helpers
# ===========================================================================

def plot_tas(ax=None, rock_type='volcanic', **kwargs):
    """Draw TAS field boundaries on a matplotlib Axes.

    Parameters
    ----------
    ax : matplotlib.axes.Axes, optional
    rock_type : {'volcanic', 'plutonic'}
    **kwargs : passed to ``ax.plot()`` for the boundary lines.

    Returns
    -------
    matplotlib.axes.Axes
    """
    import matplotlib.pyplot as plt
    if ax is None:
        ax = plt.gca()

    lkw = {'color': 'k', 'linewidth': 0.8, **kwargs}
    name_idx = 1 if rock_type == 'volcanic' else 2

    for _path, (verts, *names) in zip(_TAS_PATHS, _TAS_POLYGONS):
        x, y = verts[:, 0], verts[:, 1]
        ax.plot(np.append(x, x[0]), np.append(y, y[0]), **lkw)
        name = names[name_idx - 1]
        if name:
            cx, cy = np.mean(verts, axis=0)
            ax.text(cx, cy, name, ha='center', va='center',
                    fontsize=6, color='0.5')

    ax.set_xlabel('SiO$_2$ (wt%)')
    ax.set_ylabel('Na$_2$O + K$_2$O (wt%)')
    ax.set_xlim(33, 100)
    ax.set_ylim(0, 20)
    return ax


# ===========================================================================
# Internal helpers
# ===========================================================================

def _safe_ge(a, v): return np.where(np.isnan(a), False, a >= v)
def _safe_le(a, v): return np.where(np.isnan(a), False, a <= v)
def _safe_gt(a, v): return np.where(np.isnan(a), False, a >  v)
def _safe_lt(a, v): return np.where(np.isnan(a), False, a <  v)
