"""
Age data cleaning, resolution, and correction.

``age_from_name(name)``
    Look up a geologic time interval by name in the ICS time scale CSV,
    returning (age_min, age_mid, age_max) in Ma.

``resolve_ages(data)``
    Convert age columns that may contain numeric values *or* geologic time
    scale names (e.g. ``'Archean'``, ``'Holocene'``) into numerical Ma values,
    then run ``age_correction`` to sanitise the result.

``age_correction(data, age_var=200)``
    Sanitise age_min / age / age_max columns: fix negative values, correct
    inverted min/max, fill missing ages from min/max average, and NaN ages
    whose min–max range exceeds *age_var* Ma.
"""

from __future__ import annotations

import functools
import pathlib

import numpy as np
import pandas as pd


_DATA_DIR = pathlib.Path(__file__).resolve().parent.parent.parent / 'data'
_GTS_PATH = _DATA_DIR / 'geologic_time_scale.csv'


@functools.lru_cache(maxsize=1)
def _load_timescale() -> pd.DataFrame:
    """Load the ICS geologic time scale CSV (cached after first call)."""
    ts = pd.read_csv(_GTS_PATH)
    ts['name_lower'] = ts['name'].str.strip().str.lower()
    return ts


def age_from_name(name: str) -> tuple[float, float, float]:
    """Return (age_min, age_mid, age_max) in Ma for a geologic time interval.

    Looks up *name* (case-insensitive) in ``src/data/geologic_time_scale.csv``.
    Returns ``(NaN, NaN, NaN)`` if the name is not found.

    Parameters
    ----------
    name : str
        Geologic interval name, e.g. ``'Archean'``, ``'Holocene'``,
        ``'Cretaceous'``.

    Returns
    -------
    age_min, age_mid, age_max : float
        Ages in Ma (million years).  ``age_mid = (age_min + age_max) / 2``.

    Examples
    --------
    >>> age_from_name('Archean')
    (2500.0, 3250.0, 4000.0)
    >>> age_from_name('Holocene')
    (0.0, 0.00585, 0.0117)
    """
    ts = _load_timescale()
    key = name.strip().lower()
    row = ts[ts['name_lower'] == key]
    if row.empty:
        return np.nan, np.nan, np.nan
    age_min = float(row.iloc[0]['date_min'])
    age_max = float(row.iloc[0]['date_max'])
    return age_min, 0.5 * (age_min + age_max), age_max


def _parse_age_value(v) -> tuple[float, float, float]:
    """Convert a single age value to (min, mid, max) in Ma.

    Handles:
    - numeric (float/int)  → (NaN, value, NaN)
    - numeric string       → (NaN, value, NaN)
    - time scale name      → from age_from_name()
    - None / NaN           → (NaN, NaN, NaN)
    """
    if v is None:
        return np.nan, np.nan, np.nan
    if isinstance(v, (int, float, np.floating)):
        if np.isnan(v):
            return np.nan, np.nan, np.nan
        return np.nan, float(v), np.nan
    sv = str(v).strip()
    if not sv:
        return np.nan, np.nan, np.nan
    try:
        return np.nan, float(sv), np.nan
    except ValueError:
        return age_from_name(sv)


def resolve_ages(data: pd.DataFrame, age_var: float = 200) -> pd.DataFrame:
    """Resolve age columns containing numeric values or time scale names.

    Converts ``age``, ``age_min``, and ``age_max`` columns (which may hold
    numeric values *or* strings such as ``'Archean'`` or ``'500'``) into
    purely numerical Ma values, then runs :func:`age_correction` to sanitise
    bounds and compute ``avg_age``.

    Parameters
    ----------
    data : pandas.DataFrame
        Must contain at least one of ``'age'``, ``'age_min'``, ``'age_max'``.
    age_var : float
        Forwarded to :func:`age_correction` (default 200 Ma).

    Returns
    -------
    pandas.DataFrame
        Copy of *data* with numeric ``age``, ``age_min``, ``age_max``, and
        ``avg_age`` columns.
    """
    out = data.copy()

    for col in ('age', 'age_min', 'age_max'):
        if col not in out.columns:
            continue
        series = out[col]
        # If already all-numeric, nothing to do
        if pd.api.types.is_numeric_dtype(series):
            continue
        # Mixed or object dtype — resolve element by element
        resolved = np.full(len(series), np.nan)
        for i, v in enumerate(series):
            _, mid, _ = _parse_age_value(v)
            resolved[i] = mid
        out[col] = resolved

    return age_correction(out, age_var=age_var)


def age_correction(data, age_var=200):
    """Clean and reconcile age columns in a geochemical database.

    Applies the following fixes in order:

    1. Take absolute value of negative ages (sign errors in reporting).
    2. Convert ages > 4600 Ma from years to Ma (÷ 1 000 000).
    3. Swap inverted age_min / age_max pairs.
    4. Where age is NaN but both age_min and age_max exist, set
       age = (age_min + age_max) / 2.
    5. Fix age < age_min or age > age_max by swapping with the offending bound.
    6. Set avg_age to NaN where (age_max − age_min) > age_var.

    Parameters
    ----------
    data : pandas.DataFrame
        Must contain at least an 'age' column.  'age_min' and 'age_max' are
        used when present.
    age_var : float
        Maximum permitted range (age_max − age_min) in Ma.  Rows with a
        larger range have their avg_age set to NaN (default 200).

    Returns
    -------
    pandas.DataFrame
        Copy of *data* with an added 'avg_age' column and corrected
        age_min / age / age_max columns.
    """
    out = data.copy()

    def _col(name):
        if name in out.columns:
            return out[name].to_numpy(float)
        return np.full(len(out), np.nan)

    age     = _col('age')
    age_min = _col('age_min')
    age_max = _col('age_max')

    # 1. Absolute value — negative ages are usually sign errors
    age     = np.abs(age)
    age_min = np.abs(age_min)
    age_max = np.abs(age_max)

    # 2. Ages reported in years instead of Ma
    for arr, label in [(age, 'age'), (age_min, 'age_min'), (age_max, 'age_max')]:
        mask = arr > 4600
        n = mask.sum()
        if n:
            arr[mask] /= 1_000_000
            print(f'  age_correction: {n} {label} values > 4600 converted from years to Ma')

    # 3. Swap inverted min / max
    inv = np.isfinite(age_min) & np.isfinite(age_max) & (age_max < age_min)
    if inv.any():
        age_min[inv], age_max[inv] = age_max[inv].copy(), age_min[inv].copy()
        print(f'  age_correction: {inv.sum():,} inverted age_min/age_max pairs swapped')

    # 4. Fill missing age from average of min and max
    avg_age = age.copy()
    fill = np.isnan(avg_age) & np.isfinite(age_min) & np.isfinite(age_max)
    if fill.any():
        avg_age[fill] = (age_min[fill] + age_max[fill]) / 2
        age[fill]     = avg_age[fill]

    # 5a. Fix age < age_min
    lt_min = np.isfinite(age) & np.isfinite(age_min) & (age < age_min)
    if lt_min.any():
        age_min[lt_min], age[lt_min] = age[lt_min].copy(), age_min[lt_min].copy()
        avg_age[lt_min] = age[lt_min]
        print(f'  age_correction: {lt_min.sum():,} age < age_min corrected')

    # 5b. Fix age > age_max
    gt_max = np.isfinite(age) & np.isfinite(age_max) & (age > age_max)
    if gt_max.any():
        age_max[gt_max], age[gt_max] = age[gt_max].copy(), age_max[gt_max].copy()
        avg_age[gt_max] = age[gt_max]
        print(f'  age_correction: {gt_max.sum():,} age > age_max corrected')

    # 6. NaN avg_age where range is too wide
    wide = (np.isfinite(avg_age)
            & np.isfinite(age_min) & np.isfinite(age_max)
            & ((age_max - age_min) > age_var))
    if wide.any():
        avg_age[wide] = np.nan
        print(f'  age_correction: {wide.sum():,} ages set to NaN (range > {age_var} Ma)')

    print(f'  age_correction: {np.sum(np.isfinite(avg_age)):,} samples with valid avg_age')

    # Write corrected columns back
    if 'age'     in out.columns: out['age']     = age
    if 'age_min' in out.columns: out['age_min'] = age_min
    if 'age_max' in out.columns: out['age_max'] = age_max
    out['avg_age'] = avg_age

    return out
