"""
Classification of metamorphic rocks by facies and texture.

Uses text matching against rock name and facies description fields.
No geochemical discrimination — metamorphic grade is inferred from
reported mineralogy, not oxide composition.
"""

import numpy as np
import pandas as pd


# Metamorphic facies keywords: (facies_label, [keywords_to_match])
_FACIES = [
    ('zeolite',       ['zeolite']),
    ('hornfels',      ['hornfels', 'hornfel', 'granofel']),
    ('sanidinite',    ['sanidinite']),
    ('prehnite-pumpellyite', ['phrenite', 'pumpellyite', 'phyllit']),
    ('blueschist',    ['blueschist']),
    ('greenschist',   ['greenschist']),
    ('amphibolite',   ['amphibolite', 'amphibolit', 'hornblendite']),
    ('granulite',     ['granulite', 'granulit', 'charnock',
                       'enderbite', 'mangerite', 'jotunite',
                       'farsundite', 'opdalite']),
    ('eclogite',      ['eclogite', 'eclogit']),
]

# Metamorphic texture keywords: (texture_label, [keywords])
_TEXTURES = [
    ('slate',      ['slate', 'slaty']),
    ('schist',     ['schist']),
    ('gneiss',     ['gneiss']),
    ('migmatite',  ['migmatite', 'migmatit']),
]


def _match_keywords(text_series, keywords):
    """Return boolean array: True where any keyword is a substring of the text."""
    mask = np.zeros(len(text_series), dtype=bool)
    lower = text_series.str.lower().fillna('')
    for kw in keywords:
        mask |= lower.str.contains(kw, regex=False)
    return mask


def metamorphic_class(data, name_col='rock_name', facies_col='rock_facies'):
    """Assign metamorphic facies and texture from descriptive text fields.

    Parameters
    ----------
    data : pandas.DataFrame
    name_col : str
        Column containing the rock name (string).
    facies_col : str
        Column containing the metamorphic facies description (string).

    Returns
    -------
    pandas.DataFrame
        Input data with two new columns added:
        - 'met_facies': metamorphic facies name (empty string if unassigned)
        - 'met_texture': metamorphic texture name (empty string if unassigned)
    """
    n = len(data)
    out = data.copy()

    name   = data[name_col].astype(str)   if name_col   in data.columns else pd.Series([''] * n)
    facies = data[facies_col].astype(str) if facies_col in data.columns else pd.Series([''] * n)

    met_facies  = np.full(n, '', dtype=object)
    met_texture = np.full(n, '', dtype=object)

    for label, keywords in _FACIES:
        mask = _match_keywords(name, keywords) | _match_keywords(facies, keywords)
        met_facies[mask] = label

    for label, keywords in _TEXTURES:
        mask = _match_keywords(name, keywords) | _match_keywords(facies, keywords)
        met_texture[mask] = label

    out['met_facies']  = met_facies
    out['met_texture'] = met_texture
    return out


# ---------------------------------------------------------------------------
# Rock origin assignment for unlabelled metamorphic samples
# ---------------------------------------------------------------------------

# Keyword lists for protolith inference from rock_name
_ORIGIN_KEYWORDS = {
    'metaigneous': [
        'ingeous', 'orthogneiss', 'orthoamphibolite',   # 'ingeous' is a DB typo for 'igneous'
    ],
    'metaplutonic': [
        'granit', 'grano', 'diorit', 'tonolit', 'ttg', 'gabbro',
        'monzonit', 'syenit', 'plutonic',
        'charonockit', 'enderbite', 'norit', 'mangerite',
        'jotunite', 'farsundite', 'opdalite',
        'dike', 'dyke', 'sill',
    ],
    'metavolcanic': [
        'basalt', 'basanit', 'rhyolit', 'dacit', 'andesit', 'trachy',
        'volcanic', 'lava', 'tuff',
    ],
    'metasedimentary': [
        'sediment', 'conglomerat', 'carbon', 'graphit', 'slate', 'shale',
        'mud', 'silt', 'sand', 'lime', 'marble', 'dolomit', 'para',
    ],
}

# Priority order: metaigneous first, then more specific sub-types
_ORIGIN_ORDER = ['metaigneous', 'metaplutonic', 'metavolcanic', 'metasedimentary']


def adjust_origin(data, name_col='rock_name', group_col='rock_group',
                  origin_col='rock_origin'):
    """Infer protolith origin for unlabelled metamorphic rocks.

    Searches *name_col* for keywords that indicate an igneous or sedimentary
    protolith and updates *origin_col* accordingly.  Only samples whose
    *group_col* is 'metamorphic' and whose *origin_col* does not already have
    a specific sub-classification are updated.

    Parameters
    ----------
    data : pandas.DataFrame
    name_col : str
        Column with the rock name string (default ``'rock_name'``).
    group_col : str
        Column with the broad rock group (default ``'rock_group'``).
    origin_col : str
        Column to update with the inferred origin (default ``'rock_origin'``).

    Returns
    -------
    pandas.DataFrame
        Copy of *data* with *origin_col* updated.
    """
    out    = data.copy()
    n      = len(out)
    names  = out[name_col].astype(str).str.lower() if name_col in out.columns else pd.Series([''] * n)
    group  = out[group_col].astype(str).str.lower() if group_col in out.columns else pd.Series([''] * n)
    origin = out[origin_col].astype(str).str.lower() if origin_col in out.columns else pd.Series([''] * n)

    is_meta = group == 'metamorphic'

    for new_origin in _ORIGIN_ORDER:
        # Only update samples that are metamorphic and not already tagged as this sub-type
        candidates = is_meta & (origin != new_origin)

        # For metaplutonic / metavolcanic / metasedimentary: also exclude
        # samples already assigned 'metaigneous' (more general tag)
        if new_origin in ('metaplutonic', 'metavolcanic'):
            candidates = candidates & (origin != 'metaigneous')
        if new_origin == 'metasedimentary':
            candidates = candidates & ~origin.str.startswith('meta')

        keywords = _ORIGIN_KEYWORDS[new_origin]
        matched  = np.zeros(n, dtype=bool)
        for kw in keywords:
            matched |= names.str.contains(kw, regex=False)

        update = candidates & matched
        if update.any():
            out.loc[update, origin_col] = new_origin
            print(f'  adjust_origin: {update.sum():,} → {new_origin}')

    return out
