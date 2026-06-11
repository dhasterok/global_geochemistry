"""
Loader for the global geochemistry database.

The primary database is a large CSV exported from a PostgreSQL database.
The default path is recorded here but can be overridden at call time or
via the ``GEOCHEM_DB`` environment variable.
"""

from __future__ import annotations

import os
import pathlib
import warnings

import pandas as pd


_HARDCODED_PATH = pathlib.Path(
    os.path.expanduser(
        '~/Library/CloudStorage/GoogleDrive-dhasterok@gmail.com/'
        'My Drive/heat_production/database/export/'
        '2024_01_23/database_2024_01_23.csv'
    )
)


def _resolve_db_path() -> pathlib.Path:
    """Return database path from config/paths.yml, with hardcoded fallback."""
    try:
        import yaml
        here = pathlib.Path(__file__).resolve()
        for parent in here.parents:
            cfg = parent / 'config' / 'paths.yml'
            if cfg.exists():
                with cfg.open() as f:
                    conf = yaml.safe_load(f)
                version = conf['database']['default_version']
                raw = conf['database']['path'].format(version=version)
                return pathlib.Path(os.path.expanduser(raw))
    except Exception:
        pass
    return _HARDCODED_PATH


_DEFAULT_PATH = _resolve_db_path()


def load_database(
    path: str | pathlib.Path | None = None,
    columns: list[str] | None = None,
    rock_category: str | list[str] | None = None,
    low_memory: bool = False,
    **read_csv_kwargs,
) -> pd.DataFrame:
    """Load the global geochemistry database.

    Parameters
    ----------
    path : str or path-like, optional
        Path to the database CSV.  Falls back to the ``GEOCHEM_DB``
        environment variable, then to the default Google Drive path.
    columns : list of str, optional
        Subset of columns to load (passed as *usecols* to
        :func:`pandas.read_csv`).  Loads all columns when None.
    rock_category : str or list of str, optional
        Filter to one or more rock categories after loading
        (e.g. ``'igneous'``, ``['igneous', 'metamorphic']``).
        Requires a ``'rock_category'`` column in the database.
    low_memory : bool
        Passed to :func:`pandas.read_csv`.  Set to ``True`` for very
        large files when memory is tight (slower, may change dtypes).
    **read_csv_kwargs
        Additional keyword arguments forwarded to :func:`pandas.read_csv`.

    Returns
    -------
    pandas.DataFrame
    """
    if path is None:
        env = os.environ.get('GEOCHEM_DB')
        path = pathlib.Path(env) if env else _DEFAULT_PATH

    path = pathlib.Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"Database not found at {path}.\n"
            "Set the GEOCHEM_DB environment variable or pass path= explicitly."
        )

    kw = dict(low_memory=low_memory)
    if columns is not None:
        kw['usecols'] = columns
    kw.update(read_csv_kwargs)

    df = pd.read_csv(path, **kw)

    if rock_category is not None:
        if 'rock_category' not in df.columns:
            warnings.warn("'rock_category' column not found; ignoring filter.")
        else:
            if isinstance(rock_category, str):
                rock_category = [rock_category]
            mask = df['rock_category'].str.lower().isin(
                [r.lower() for r in rock_category]
            )
            df = df[mask].reset_index(drop=True)

    return df
