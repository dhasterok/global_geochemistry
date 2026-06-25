"""
EarthChem portal data reader.

Reads the Excel or CSV files produced by earthchem.org and returns a
tidy :class:`pandas.DataFrame` with standardised column names matching
the project's database convention.

Column normalisation applied:
- Lower-case, spaces → underscores, subscript Unicode → ASCII digits
- Alias mapping (e.g. ``"sample id"`` → ``"sample_name"``)
- BDL tokens (``"bdl"``, ``"n.d."``, …) → 0.0
- ``"<value"`` → negative (BDL convention), ``">"`` → strip prefix
- Unit-suffix stripping with conversion (``_ppb`` × 10³, ``_ppt`` × 10⁻⁶)
- Elemental ppm columns get ``_ppm`` suffix normalised
"""

from __future__ import annotations

import pathlib
import re

import numpy as np
import pandas as pd

from global_geochemistry.geochem.fileio.pgcsv import (
    _clean_name,
    _resolve_name,
    _apply_units,
    _BDL_TOKENS,
    _ELEMENTS,
)


def ecread(filename) -> pd.DataFrame:
    """Read an EarthChem portal export file.

    Accepts ``.xls``, ``.xlsx``, or ``.csv``.  Column names are normalised
    to the project database convention (lower-case, underscores, aliases
    resolved).  Numeric columns containing BDL tokens or ``<value`` notation
    are converted using the same conventions as :func:`~.pgcsv.pgcsv`.

    Parameters
    ----------
    filename : str or path-like

    Returns
    -------
    pandas.DataFrame
    """
    path = pathlib.Path(filename)
    suffix = path.suffix.lower()

    if suffix in ('.xls', '.xlsx'):
        df = pd.read_excel(path, header=0, dtype=object)
    elif suffix == '.csv':
        df = pd.read_csv(path, header=0, low_memory=False, dtype=object)
    else:
        raise ValueError(f"Unsupported file type: '{suffix}'")

    df = _normalise_columns(df)
    df = _convert_values(df)
    return df


def _normalise_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Clean and alias column names, applying unit conversion where encoded."""
    rename_map: dict[str, str] = {}
    unit_factors: dict[str, float] = {}

    for raw in df.columns:
        name = _clean_name(str(raw))
        name = _resolve_name(name)

        if name.endswith('_ppb'):
            unit_factors[raw] = 1e3
            name = name[:-4]
        elif name.endswith('_ppt'):
            unit_factors[raw] = 1e-6
            name = name[:-4]
        elif name.endswith('_wtper'):
            unit_factors[raw] = 1e4
            name = name[:-6]
        elif name.endswith('_per'):
            unit_factors[raw] = 1e4
            name = name[:-4]
        elif name.endswith('_ppm'):
            name = name[:-4]

        bare = re.sub(r'_(ppb|ppt|per|wtper|ppm)$', '', name)
        if bare in _ELEMENTS:
            name = bare + '_ppm'

        rename_map[raw] = name

    for raw_col, factor in unit_factors.items():
        df[raw_col] = pd.to_numeric(df[raw_col], errors='coerce') * factor

    df = df.rename(columns=rename_map)

    # Merge duplicate column names by averaging numeric values
    if df.columns.duplicated().any():
        processed: list[str] = []
        seen: set[str] = set()
        for col in df.columns.unique():
            if col in seen:
                continue
            seen.add(col)
            dupes = df.loc[:, df.columns == col]
            if dupes.shape[1] > 1:
                numeric = dupes.apply(pd.to_numeric, errors='coerce')
                df[col + '_merged_'] = numeric.mean(axis=1)
                processed.append(col + '_merged_')
            else:
                processed.append(col)
        df = df.loc[:, ~df.columns.duplicated(keep=False)]
        merged = {c: c[:-8] for c in df.columns if c.endswith('_merged_')}
        if merged:
            df = df.rename(columns=merged)

    return df


def _convert_values(df: pd.DataFrame) -> pd.DataFrame:
    """Apply BDL and ``<value`` conversion to numeric-looking columns."""
    for col in df.columns:
        series = df[col]
        sample = series.dropna().astype(str).str.strip().str.lower()

        is_element = col.endswith('_ppm') or col.rstrip('_ppm') in _ELEMENTS
        has_bdl = sample.isin(_BDL_TOKENS).any()
        has_lt = sample.str.startswith('<').any()
        looks_numeric = (
            is_element or has_bdl or has_lt
            or sample.str.match(r'^[-+]?\d').any()
        )
        if not looks_numeric:
            continue

        cleaned: list = []
        all_numeric = True
        for v in series:
            sv = str(v).strip().lower() if v is not None else 'nan'
            if sv in ('nan', '', 'none'):
                cleaned.append(np.nan)
            elif sv in _BDL_TOKENS:
                cleaned.append(0.0)
            elif sv == '-':
                cleaned.append(np.nan)
            elif sv.startswith('<'):
                try:
                    cleaned.append(-float(sv[1:].strip()))
                except ValueError:
                    cleaned.append(np.nan)
            elif sv.startswith('>'):
                try:
                    cleaned.append(float(sv[1:].strip()))
                except ValueError:
                    cleaned.append(np.nan)
            else:
                try:
                    cleaned.append(float(v))
                except (ValueError, TypeError):
                    cleaned.append(v)
                    all_numeric = False

        if all_numeric:
            df[col] = pd.array(cleaned, dtype=object)
            try:
                df[col] = df[col].astype(float)
            except (ValueError, TypeError):
                pass

    return df
