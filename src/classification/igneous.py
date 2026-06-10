"""
Classification of igneous rocks.

Implements:
- TAS (Total Alkali–Silica) classification for volcanic and plutonic rocks
  based on Le Bas et al. (1986) and Middlemost (1994).
- High-Mg volcanic and ultramafic rock types (Le Bas & Streckeisen 1991).
- Mg number computation.

References:
    Le Bas, M.J., Le Maitre, R.W., Streckeisen, A., Zanettin, B. (1986)
    A chemical classification of volcanic rocks based on the total alkali
    silica diagram.  Journal of Petrology 27(3), 745–750.

    Middlemost, E.A.K. (1994) Naming materials in the magma/igneous rock
    system.  Earth-Science Reviews 37, 215–224.

    Le Bas, M.J. & Streckeisen, A.L. (1991) The IUGS systematics of igneous
    rocks.  Journal of the Geological Society London 148, 825–833.
"""

import numpy as np
import pandas as pd
from matplotlib.path import Path


# ---------------------------------------------------------------------------
# TAS polygon definitions (SiO2, Na2O+K2O) — from Middlemost (1994)
# Each entry: (vertices_array, volcanic_name, plutonic_name)
# ---------------------------------------------------------------------------

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

# Pre-build matplotlib Path objects for fast point-in-polygon testing
_TAS_PATHS = [Path(v) for v, _, _ in _TAS_POLYGONS]


def tas(data, rock_category=None):
    """Classify igneous rocks using the Total Alkali–Silica diagram.

    Parameters
    ----------
    data : pandas.DataFrame
        Must contain 'sio2', 'na2o', 'k2o', 'mgo', 'feo_tot', 'tio2'
        columns in wt%.
    rock_category : array-like of str, optional
        Per-sample classification as 'volcanic', 'plutonic', or other/NaN.
        If None, all samples are classified with volcanic names.  Samples
        with category 'plutonic' receive plutonic names.

    Returns
    -------
    pandas.Series
        Rock names (string), empty string where classification could not
        be assigned.

    Notes
    -----
    High-Mg rocks (MgO ≥ 12 wt%) receive overriding names.
    Ultramafic plutonic rocks are further subdivided by TiO2 content.
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
    pts_valid = np.column_stack([sio2[valid], ta[valid]])

    # Determine rock category masks
    is_volcanic = np.zeros(n, dtype=bool)
    is_plutonic = np.zeros(n, dtype=bool)
    if rock_category is None:
        is_volcanic = valid.copy()  # classify all with volcanic names when unspecified
    else:
        rc = np.asarray(rock_category, dtype=object)
        is_volcanic = np.isin(rc, ['volcanic'])
        is_plutonic = np.isin(rc, ['plutonic'])

    for path, (_, volc_name, plut_name) in zip(_TAS_PATHS, _TAS_POLYGONS):
        inside = np.zeros(n, dtype=bool)
        if pts_valid.shape[0]:
            inside[valid] = path.contains_points(pts_valid)
        names[inside & is_volcanic] = volc_name
        names[inside & is_plutonic] = plut_name

    # --- High-Mg overrides ---
    # Applied to both volcanic and plutonic unless otherwise noted.
    mafic_field = valid & (sio2 >= 33) & (sio2 < 65)

    # Boninite/sanukitoid
    mask = mafic_field & (sio2 >= 52) & _safe_ge(mgo, 8) & _safe_lt(tio2, 0.5)
    names[mask & is_volcanic] = 'boninite'
    names[mask & is_plutonic] = 'sanukitoid'

    # Picrite / ferropicrite variants — no volcanic/plutonic distinction in MATLAB
    mask = mafic_field & (sio2 < 52) & _safe_ge(mgo, 12) & (ta < 3)
    names[mask & _safe_ge(feo, 13)] = 'ferropicrite'
    names[mask & _safe_lt(feo, 13)] = 'picrite'

    mask = mafic_field & (sio2 < 52) & _safe_ge(mgo, 12) & (ta >= 3)
    names[mask & _safe_ge(feo, 13)] = 'alkali ferropicrite'
    names[mask & _safe_lt(feo, 13)] = 'alkali picrite'

    # Meimechite (very high MgO, high TiO2)
    mask = mafic_field & (sio2 < 52) & _safe_ge(mgo, 18) & _safe_gt(tio2, 1) & (ta < 2)
    names[mask] = 'meimechite'

    # Komatiite
    mask = mafic_field & (sio2 < 45) & _safe_ge(mgo, 18) & _safe_lt(tio2, 1) & (ta < 2)
    names[mask & is_volcanic] = 'komatiite'
    names[mask & is_plutonic] = 'intrusive komatiite'

    # Basaltic komatiite
    mask = mafic_field & (sio2 >= 45) & (sio2 < 52) & _safe_ge(mgo, 18) & _safe_lt(tio2, 1) & (ta < 2)
    names[mask & is_volcanic] = 'basaltic komatiite'
    names[mask & is_plutonic] = 'gabbroic komatiite'

    # Plutonic ultramafic subdivision by TiO2 (Reverdatto et al. 2008)
    mask = is_plutonic & valid & (sio2 >= 33) & (sio2 < 45) & _safe_gt(tio2, 0) & _safe_le(tio2, 0.3)
    names[mask & _safe_ge(mgo, 35)]                          = 'mantle peridotite'
    names[mask & _safe_ge(mgo, 18) & _safe_lt(mgo, 35)]     = 'mantle pyroxenite'

    return pd.Series(names, index=data.index, dtype=object)


def mgnum(mgo, feo_tot, fe3_fraction=0.0):
    """Compute Mg number: Mg/(Mg + Fe²⁺) molar.

    Parameters
    ----------
    mgo : array-like
        MgO concentration in wt%.
    feo_tot : array-like
        Total iron as FeO in wt%.
    fe3_fraction : float or array-like
        Fraction of total iron as Fe³⁺ (0–1).  Default 0 (all Fe as Fe²⁺).

    Returns
    -------
    numpy.ndarray
    """
    from src.utils.molecular import MolecularWeightCalculator
    _mwc = MolecularWeightCalculator()

    mgo     = np.asarray(mgo,     dtype=float)
    feo_tot = np.asarray(feo_tot, dtype=float)
    fe3_frac = np.asarray(fe3_fraction, dtype=float)

    n_mg  = mgo / _mwc.molecular_weight('MgO')
    n_feo = feo_tot * (1 - fe3_frac) / _mwc.molecular_weight('FeO')
    with np.errstate(invalid='ignore', divide='ignore'):
        mg = n_mg / (n_mg + n_feo)
    mg = np.where(np.isfinite(mg), mg, np.nan)
    return mg


def plot_tas(ax=None, rock_type='volcanic', **kwargs):
    """Draw TAS field boundaries on a matplotlib axes.

    Parameters
    ----------
    ax : matplotlib.axes.Axes, optional
        Target axes; uses current axes if None.
    rock_type : {'volcanic', 'plutonic'}
        Which rock names to annotate.
    **kwargs
        Passed to ``ax.plot()`` for boundary lines.

    Returns
    -------
    matplotlib.axes.Axes
    """
    import matplotlib.pyplot as plt

    if ax is None:
        ax = plt.gca()

    lkw = {'color': 'k', 'linewidth': 0.8, **kwargs}
    name_idx = 1 if rock_type == 'volcanic' else 2

    for path, (verts, *names) in zip(_TAS_PATHS, _TAS_POLYGONS):
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


# ---------------------------------------------------------------------------
# Helper: NaN-safe comparisons
# ---------------------------------------------------------------------------

def _safe_ge(a, v): return np.where(np.isnan(a), False, a >= v)
def _safe_le(a, v): return np.where(np.isnan(a), False, a <= v)
def _safe_gt(a, v): return np.where(np.isnan(a), False, a >  v)
def _safe_lt(a, v): return np.where(np.isnan(a), False, a <  v)
