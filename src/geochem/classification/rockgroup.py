"""
Rock-group filtering.

Generalizes the rock-type/rock-origin boolean queries used internally
throughout this package (e.g. the private igneous/sedimentary/metamorphic
masks in `geochem.fileio.pipeline`) into a single, publicly reusable
function covering the full type vocabulary of MATLAB's `rockgroup.m`.
"""

import warnings

import pandas as pd


def _eq(data, column, *values):
    if column not in data.columns:
        return pd.Series(False, index=data.index)
    lowered = data[column].astype(str).str.lower()
    return lowered.isin([v.lower() for v in values])


def _igneous_protolith(data):
    base = _eq(data, 'rock_group', 'igneous') | _eq(data, 'rock_origin', 'metaigneous', 'metaplutonic', 'metavolcanic')
    if 'protolith_est' not in data.columns:
        warnings.warn('Asked for estimated protoliths, but they have not yet been determined.')
        return base
    return base | (
        _eq(data, 'rock_group', 'metamorphic')
        & _eq(data, 'protolith_est', 'metaigneous')
        & ~_eq(data, 'rock_origin', 'metasedimentary')
    )


def _sedimentary_protolith(data):
    base = _eq(data, 'rock_group', 'sedimentary') | _eq(data, 'rock_origin', 'metasedimentary')
    if 'protolith_est' not in data.columns:
        warnings.warn('Asked for estimated protoliths, but they have not yet been determined.')
        return base
    return base | (
        _eq(data, 'rock_group', 'metamorphic')
        & _eq(data, 'protolith_est', 'metasedimentary')
        & ~_eq(data, 'rock_origin', 'metaigneous', 'metaplutonic', 'metavolcanic')
    )


_TYPES = {
    'igneous':               lambda d: _eq(d, 'rock_group', 'igneous'),
    'all igneous':           lambda d: _eq(d, 'rock_group', 'igneous')
                                        | _eq(d, 'rock_origin', 'metaigneous', 'metaplutonic', 'metavolcanic'),
    'igneous protolith':     _igneous_protolith,
    'volcanic':              lambda d: _eq(d, 'rock_origin', 'volcanic'),
    'metavolcanic':          lambda d: _eq(d, 'rock_origin', 'metavolcanic'),
    'all volcanic':          lambda d: _eq(d, 'rock_origin', 'volcanic', 'metavolcanic'),
    'plutonic':              lambda d: _eq(d, 'rock_origin', 'plutonic'),
    'metaplutonic':          lambda d: _eq(d, 'rock_origin', 'metaplutonic'),
    'all plutonic':          lambda d: _eq(d, 'rock_origin', 'plutonic', 'metaplutonic'),
    'sedimentary':           lambda d: _eq(d, 'rock_group', 'sedimentary'),
    'all seds':              lambda d: _eq(d, 'rock_group', 'sedimentary') | _eq(d, 'rock_origin', 'metasedimentary'),
    'sedimentary protolith': _sedimentary_protolith,
    'metamorphic':           lambda d: _eq(d, 'rock_group', 'metamorphic'),
    'metaigneous':           lambda d: _eq(d, 'rock_origin', 'metaigneous', 'metaplutonic', 'metavolcanic'),
    'metasedimentary':       lambda d: _eq(d, 'rock_origin', 'metasedimentary'),
    'mineral':               lambda d: _eq(d, 'material', 'mineral'),
}


def rockgroup(data, rock_type):
    """Returns a boolean mask selecting samples of the desired rock group.

    Parameters
    ----------
    data : pandas.DataFrame
        Must contain some subset of ``rock_group``, ``rock_origin``,
        ``protolith_est``, ``material``, depending on `rock_type`.
    rock_type : str
        One of: ``'igneous'``, ``'all igneous'``, ``'igneous protolith'``,
        ``'volcanic'``, ``'metavolcanic'``, ``'all volcanic'``,
        ``'plutonic'``, ``'metaplutonic'``, ``'all plutonic'``,
        ``'sedimentary'``, ``'all seds'``, ``'sedimentary protolith'``,
        ``'metamorphic'``, ``'metaigneous'``, ``'metasedimentary'``,
        ``'mineral'``.

    Returns
    -------
    pandas.Series of bool
    """
    key = rock_type.lower()
    if key not in _TYPES:
        raise ValueError(f"Unknown rock_type: {rock_type!r}. Choose from {sorted(_TYPES)}.")
    return _TYPES[key](data)
