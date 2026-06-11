"""
EarthChem data reader.

Reads the Excel or CSV files produced by earthchem.org and returns a
tidy :class:`pandas.DataFrame` with lower-case column names.
"""

import pathlib
import pandas as pd


def ecread(filename):
    """Read an EarthChem export file.

    Accepts ``.xls``, ``.xlsx``, or ``.csv``.  Column names are
    lower-cased and stripped of leading/trailing whitespace.

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
        df = pd.read_excel(path, header=0)
    elif suffix == '.csv':
        df = pd.read_csv(path, header=0, low_memory=False)
    else:
        raise ValueError(f"Unsupported file type: '{suffix}'")

    df.columns = [str(c).strip().lower() for c in df.columns]
    return df
