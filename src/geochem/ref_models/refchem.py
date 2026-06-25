"""
Reference geochemistry reservoir loader.

Reads ``earthref.xlsx`` from ``src/data/`` and optionally applies the
standard processing pipeline (Fe correction, geochemical indices, physical
property estimates).

Main functions
--------------
``load_earthref(physprop=True)``
    Full pipeline — equivalent to MATLAB ``loaderef.m``.

``load_refchem(models=None, fefix=False, physprop=False)``
    Filtered load — equivalent to MATLAB ``loadrefchem.m``.
"""

from __future__ import annotations

import pathlib
import warnings

import numpy as np
import pandas as pd

_DATA_DIR = pathlib.Path(__file__).parent.parent.parent / 'data'
_EARTHREF_PATH = _DATA_DIR / 'earthref.xlsx'


def _iron_fix(df: pd.DataFrame) -> pd.DataFrame:
    """Apply Fe correction (wrapper around src.processing.iron.fefix)."""
    try:
        from global_geochemistry.geochem.processing.iron import fefix
        return fefix(df)
    except ImportError:
        warnings.warn("src.processing.iron not found; skipping Fe fix.")
        return df


def _geochem_index(df: pd.DataFrame) -> pd.DataFrame:
    try:
        from global_geochemistry.geochem.classification.indices import geochem_index
        return geochem_index(df)
    except ImportError:
        warnings.warn("src.classification.indices not found; skipping geochem indices.")
        return df


def _vpest(df: pd.DataFrame) -> pd.DataFrame:
    try:
        from global_geochemistry.geochem.physprop.seismic import vpest
        return vpest(df)
    except ImportError:
        warnings.warn("src.physprop.seismic not found; skipping Vp estimate.")
        return df


def _densest(df: pd.DataFrame) -> pd.DataFrame:
    try:
        from global_geochemistry.geochem.physprop.density import densest
        return densest(df)
    except ImportError:
        warnings.warn("src.physprop.density not found; skipping density estimate.")
        return df


def _hpest(df: pd.DataFrame) -> pd.DataFrame:
    try:
        from global_geochemistry.geochem.physprop.heat import hpest
        return hpest(df)
    except ImportError:
        warnings.warn("src.physprop.heat not found; skipping HP estimate.")
        return df


def _tcest(df: pd.DataFrame) -> pd.DataFrame:
    try:
        from global_geochemistry.geochem.physprop.thermal import tcest
        return tcest(df)
    except ImportError:
        warnings.warn("src.physprop.thermal not found; skipping Tc estimate.")
        return df


def _null_sigma_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Set estimated property columns to NaN where sigma flag is set."""
    if 'sigma' not in df.columns:
        return df
    sigma = df['sigma'].astype(bool)
    null_cols = [
        'p_velocity', 's_velocity',
        'heat_production_mass', 'heat_production',
        'density_bk', 'thermal_conductivity',
    ]
    for col in null_cols:
        if col in df.columns:
            df.loc[sigma, col] = np.nan
    return df


def load_earthref(
    path: str | pathlib.Path | None = None,
    physprop: bool = True,
) -> pd.DataFrame:
    """Load and process the earthref.xlsx reference reservoir table.

    Equivalent to MATLAB ``loaderef.m``.

    Parameters
    ----------
    path : str or path-like, optional
        Path to ``earthref.xlsx``.  Defaults to ``src/data/earthref.xlsx``.
    physprop : bool
        When True (default) compute seismic velocity, density, heat
        production, and thermal conductivity for each reservoir.

    Returns
    -------
    pandas.DataFrame
    """
    p = pathlib.Path(path) if path else _EARTHREF_PATH
    if not p.exists():
        raise FileNotFoundError(f"earthref.xlsx not found at {p}")

    df = pd.read_excel(p)
    df.columns = [str(c).strip().lower() for c in df.columns]

    df = _iron_fix(df)
    df = _geochem_index(df)

    if physprop:
        df = _vpest(df)
        df = _densest(df)
        df = _hpest(df)
        df = _tcest(df)
        df = _null_sigma_columns(df)

    return df


def load_refchem(
    models: str | list[str] | None = None,
    fefix: bool = False,
    physprop: bool = False,
    path: str | pathlib.Path | None = None,
) -> pd.DataFrame:
    """Load reference chemistry, optionally filtered to specific models.

    Equivalent to MATLAB ``loadrefchem.m``.

    Parameters
    ----------
    models : str or list of str, optional
        Filter to one or more model names (matched against the ``'model'``
        column).  Returns all rows when None.
    fefix : bool
        Apply Fe oxidation correction (default False).
    physprop : bool
        Compute physical property estimates (default False).
    path : str or path-like, optional
        Path to ``earthref.xlsx``.

    Returns
    -------
    pandas.DataFrame
    """
    p = pathlib.Path(path) if path else _EARTHREF_PATH
    if not p.exists():
        raise FileNotFoundError(f"earthref.xlsx not found at {p}")

    df = pd.read_excel(p)
    df.columns = [str(c).strip().lower() for c in df.columns]

    if models is not None and 'model' in df.columns:
        if isinstance(models, str):
            models = [models]
        models_lower = [m.lower() for m in models]
        df = df[df['model'].str.lower().isin(models_lower)].reset_index(drop=True)

    if fefix:
        df = _iron_fix(df)

    if physprop:
        df = _geochem_index(df)
        df = _vpest(df)
        df = _densest(df)
        df = _hpest(df)
        df = _null_sigma_columns(df)

    return df
