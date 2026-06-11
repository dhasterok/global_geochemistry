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
