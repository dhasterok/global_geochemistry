"""
IUGS QAP modal classification of plutonic and volcanic rocks.

Uses normative mineral proportions (CIPW) to define Q (quartz),
A (alkali feldspar), P (plagioclase), and F (feldspathoids) fractions,
then assigns rock names via polygon testing on either the Q-A-P or A-P-F
ternary face of the IUGS double triangle.

Reference:
    Le Maitre, R.W. (ed.) (2002) Igneous Rocks: A Classification and
    Glossary of Terms.  Cambridge University Press.  2nd edition.
"""

import numpy as np
import pandas as pd
from matplotlib.path import Path


# ---------------------------------------------------------------------------
# QAP polygon definitions.
# Each vertex row: [Q, A, P, F] as percent of Q+A+P+F.
# Polygons on the Q-A-P face have F = 0.
# Polygons on the A-P-F face have Q = 0.
# Source: load_qapgons.m
# ---------------------------------------------------------------------------
_QAP_POLYGONS = [
    # Q-A-P face (F = 0)  — columns: [Q, A, P, F]
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
    # Q-A-P: quartz-bearing syenites/trachytes
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
    # A-P only (Q = 0, F = 0)
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
    # A-P-F face (Q = 0)
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


def _ternary_path(verts, face):
    """Build a matplotlib Path in 2D ternary coordinates for one face.

    face : 'QAP'  → project onto (A, Q) axis ignoring P (since Q+A+P = 100)
           'APF'  → project onto (A, F) axis ignoring P
    """
    if face == 'QAP':
        # For Q-A-P: x = A/(Q+A+P), y = Q/(Q+A+P); since F=0, total = Q+A+P
        total = verts[:, 0] + verts[:, 1] + verts[:, 2]
        with np.errstate(invalid='ignore', divide='ignore'):
            x = verts[:, 1] / total   # A fraction
            y = verts[:, 0] / total   # Q fraction
    else:  # APF
        # For A-P-F: x = A/(A+P+F), y = F/(A+P+F); since Q=0, total = A+P+F
        total = verts[:, 1] + verts[:, 2] + verts[:, 3]
        with np.errstate(invalid='ignore', divide='ignore'):
            x = verts[:, 1] / total   # A fraction
            y = verts[:, 3] / total   # F fraction
    return Path(np.column_stack([x, y]))


# Determine which face each polygon belongs to and pre-build paths
_POLYGON_DATA = []
for verts, plut_name, volc_name in _QAP_POLYGONS:
    has_q = verts[:, 0].max() > 0.5   # has nonzero Q
    has_f = verts[:, 3].max() > 0.5   # has nonzero F
    face = 'APF' if has_f else 'QAP'
    _POLYGON_DATA.append({
        'face': face,
        'path': _ternary_path(verts, face),
        'plutonic': plut_name,
        'volcanic': volc_name,
    })


def qap(cipw_result, is_volcanic=None):
    """Assign QAP rock names from CIPW normative mineral proportions.

    Parameters
    ----------
    cipw_result : pandas.DataFrame
        Output of :func:`src.classification.cipw.cipwnorm`.  Must contain
        columns ``cipw_quartz``, ``cipw_orthoclase``, ``cipw_albite``,
        ``cipw_anorthite``, ``cipw_nepheline``, ``cipw_leucite``,
        ``cipw_kaliophilite``.
    is_volcanic : array-like of bool, optional
        True → volcanic name, False → plutonic.  Defaults to False (plutonic).

    Returns
    -------
    pandas.Series
        QAP rock name per sample; empty string where classification is
        not possible (Q+A+P+F = 0 or missing).
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
        Pn = np.where(valid, P / total, np.nan)
        Fn = np.where(valid, F / total, np.nan)

    # Samples with F > 10% go on the A-P-F face; others on Q-A-P
    use_apf = valid & (Fn > 0.1)
    use_qap = valid & ~use_apf

    # Build 2D test points for each face
    pts_qap = np.column_stack([An, Qn])   # (A, Q) in [0,1]
    pts_apf = np.column_stack([An, Fn])   # (A, F) in [0,1]

    names = np.full(n, '', dtype=object)

    for pd_entry in _POLYGON_DATA:
        face   = pd_entry['face']
        path   = pd_entry['path']
        p_name = pd_entry['plutonic']
        v_name = pd_entry['volcanic']

        use = use_apf if face == 'APF' else use_qap
        if not use.any():
            continue

        inside = np.zeros(n, dtype=bool)
        inside[use] = path.contains_points(
            pts_apf[use] if face == 'APF' else pts_qap[use]
        )

        names[inside & ~is_volcanic] = p_name
        names[inside &  is_volcanic] = v_name

    return pd.Series(names, index=cipw_result.index, dtype=object)
