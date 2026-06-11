"""
Convert geochemical Excel template files to database-ready CSV.

Implements the MATLAB ``pgcsv.m`` / ``parsexls`` workflow in Python.

The Excel templates have a **transposed** layout:
- Each **row** is one variable (analyte, metadata field).
- Each **column** is one sample.
- The first column contains the variable name; subsequent columns contain
  values.
- Data begins at the row whose first cell is ``'reference'``; all rows
  above it (metadata about the submission) are skipped.

Processing steps applied to every file:

1. Locate the ``'reference'`` row and use everything below it as the data
   table.
2. Normalise variable names: lower-case, whitespace → ``_``, strip
   subscript Unicode characters, fix known OCR errors and naming variants.
3. Convert units where the suffix encodes the unit:
   ``_ppb`` → ×10³, ``_ppt`` → ×10⁻⁶, ``_per`` / ``_wtper`` → ×10⁴.
4. Map below-detection-limit tokens (``bdl``, ``b.d.l.``, ``nd``, …) to
   zero; ``<value`` → negative; bare ``-`` → NaN.
5. Merge duplicate variable names by averaging numeric values.
6. Ensure required metadata columns are present (filled with empty string
   if missing).
7. Write a UTF-8 CSV to *output_dir* (one CSV per input file).

Usage
-----
Process a single file::

    pgcsv('path/to/file.xlsx', output_dir='outdir/')

Process an entire directory::

    pgcsv('path/to/spreadsheets/', output_dir='outdir/')
"""

from __future__ import annotations

import pathlib
import re
import unicodedata
import warnings

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Column name aliases (lower-case → canonical name)
# ---------------------------------------------------------------------------

_ALIASES: dict[str, str] = {
    'reference':            'bibtex',
    'bibkey':               'bibtex',
    'or':                   'orth',       # SQL reserved word
    'group':                'chem_group', # SQL reserved word
    'sfeo':    'feo_tot', 'feot':  'feo_tot', 'feo_t':  'feo_tot',
    'tfeo':    'feo_tot', 'feo*':  'feo_tot', 'feotot': 'feo_tot',
    'sfe2o3':  'fe2o3_tot', 'fe2o3t': 'fe2o3_tot',
    'tfe2o3':  'fe2o3_tot', 'fe2o3_t': 'fe2o3_tot', 'fe2o3tot': 'fe2o3_tot',
    's_tot':   's_wtper',
    'c_tot':   'c_wtper',
    'h2o+':    'h2o_plus',
    'h2o-':    'h2o_minus',
    'h2o?':    'h2o_minus',
    'geologic_province':      'prov_name',
    'subprovince_or_terrane': 'subprov1',
    'lithologic_unit':        'lith_name',
    'rock_type':              'rock_name',
    'rock type':              'rock_name',
    'rock facies':            'rock_composition',
    'rock_facies':            'rock_composition',
    'mage_ma':     'mage',
    'mage_min_ma': 'mage_min',
    'mage_max_ma': 'mage_max',
    'mage_sd_ma':  'mage_sd',
    'sum':         'total',
    'sample_id':   'sample_name',
    # OCR / typo errors
    'na20': 'na2o', '2o': 'na2o', '20': 'na2o',
    'k20':  'k2o',
    'a1203': 'al2o3', 'al203': 'al2o3',
    'si02':  'sio2',
    'ti02':  'tio2',
    'p205':  'p2o5',
    'ca0':   'cao',
    'mg0':   'mgo',
    'mn0':   'mno',
    'fe0':   'feo',
    'fe203': 'fe2o3', 'fo203': 'fe2o3',
    # Mg number aliases
    'mgo#': 'mg_num', 'mg#': 'mg_num', 'mg_#': 'mg_num',
    'mg_number': 'mg_num', 'mg number': 'mg_num',
    'mg #': 'mg_num', '#mg': 'mg_num', 'mg-no': 'mg_num',
    # LOI
    'l.o.i.': 'loi', 'l.o.i': 'loi', 'loss': 'loi',
    # Isotope ratios
    '87rb_86sr': 'rb87_sr86', '87sr_86sr': 'sr87_sr86',
    '84sr_86sr': 'sr84_sr86', '84sr_88sr': 'sr84_sr88',
    '147sm_144nd': 'sm147_nd144', '147sm_143nd': 'sm147_nd143',
    '143nd_144nd': 'nd143_nd144', '146nd_145nd': 'nd146_nd145',
    '176hf_177hf': 'hf176_hf177', '176lu_177hf': 'lu176_hf177',
    '187re_188os': 're187_os188', '187os_188os': 'os187_os188',
    '206pb_204pb': 'pb206_pb204', '207pb_204pb': 'pb207_pb204',
    '208pb_204pb': 'pb208_pb204',
    '206pb_238u': 'pb206_u238', '207pb_238u': 'pb207_u238',
    '206pb_235u': 'pb206_u235', '207pb_235u': 'pb207_u235',
    '238u_232th': 'u238_th232', '230th_232th': 'th230_th232',
    '230th_238u': 'th230_u238', '226ra_230th': 'ra226_th230',
    '238u_204pb': 'u238_pb204', '232th_204pb': 'th232_pb204',
}

# Uncertainty columns follow the same pattern with '_uncertainty' suffix
_ALIASES.update({
    k + '_uncertainty': v + '_uncertainty'
    for k, v in _ALIASES.items()
    if not k.endswith('_uncertainty')
})

# Elemental symbols (ppm columns after unit stripping)
_ELEMENTS = {
    'li','na','k','rb','cs','be','mg','ca','sr','ba','sc','y',
    'la','ce','pr','nd','sm','eu','gd','tb','dy','ho','er','tm','yb','lu',
    'ti','zr','hf','v','nb','ta','cr','mo','w','mn','fe','co','ni',
    're','ru','os','rh','pd','pt','ir','cu','ag','au',
    'zn','cd','hg','b','al','ga','in','tl','c','si','ge','sn','pb',
    'n','p','as','sb','bi','s','se','te','f','cl','br','i',
    'u','th',
}

# Numeric field names (wt%)
_NUMFIELDS = {
    'sio2','tio2','al2o3','fe2o3','fe2o3_tot','feo','feo_tot',
    'mgo','cao','na2o','k2o','p2o5','mno','sro','bao','cr2o3','nio',
    'caco3','so3','h2o_tot','h2o_plus','h2o_minus','co2','loi','total',
    'latitude','longitude','mswd','mmswd',
}

# Time fields
_TIMEFIELDS = {'age','age_min','age_max','age_sd','mage','mage_min','mage_max','mage_sd'}

# Required metadata columns (added as blank strings if missing)
_REQUIRED = {
    'author','title','year','journal','doi','bibtex',
    'lith_name','lith_alt_name','prov_name',
}

# Below-detection-limit tokens → 0
_BDL_TOKENS = {'bdl','b.d.l.','b.d.l','nd','n.d.','bd','b.d.','<dl','<d.l.'}

# Subscript Unicode → ASCII digits
_SUBSCRIPTS = str.maketrans({
    '₀':'0','₁':'1','₂':'2','₃':'3','₄':'4',
    '₅':'5','₆':'6','₇':'7','₈':'8','₉':'9',
})


def _clean_name(raw: str) -> str:
    """Normalise a raw variable name from the Excel template."""
    s = str(raw).strip()
    s = s.translate(_SUBSCRIPTS)
    for _bad in ('\r', '\x85', '\u0338', '\u2028', '\u2029'):
        s = s.replace(_bad, ' ')
    s = s.replace('/', '_').replace('\\', '_')
    paren = s.find('(')
    if paren != -1:
        s = s[:paren]
    s = re.sub(r'[^a-zA-Z0-9 _.#+*\-?<>]', '', s)
    s = s.lower().strip()
    s = re.sub(r'\s+', '_', s)
    s = s.rstrip('_')
    return s


def _resolve_name(name: str) -> str:
    """Apply aliases and strip unit suffixes from a cleaned variable name."""
    # Strip _uncertainty before alias lookup, then restore
    is_unc = name.endswith('_uncertainty')
    base = name[:-len('_uncertainty')] if is_unc else name
    base = _ALIASES.get(base, base)
    return (base + '_uncertainty') if is_unc else base


def _apply_units(name: str, series: pd.Series) -> tuple[str, pd.Series]:
    """Strip unit suffix from *name* and convert *series* accordingly.

    Returns the updated ``(name, series)`` pair.
    """
    s = series.copy()
    if name.endswith('_ppb'):
        s = pd.to_numeric(s, errors='coerce') * 1e3
        name = name[:-4]
    elif name.endswith('_ppt'):
        s = pd.to_numeric(s, errors='coerce') * 1e-6
        name = name[:-4]
    elif name.endswith('_per'):
        s = pd.to_numeric(s, errors='coerce') * 1e4
        name = name[:-4]
    elif name.endswith('_wtper'):
        s = pd.to_numeric(s, errors='coerce') * 1e4
        name = name[:-6]
    elif name.endswith('_ppm'):
        # already in ppm, just strip suffix
        s = pd.to_numeric(s, errors='coerce')
        name = name[:-4]
    return name, s


def _is_numeric_field(name: str) -> bool:
    """True if this variable should be treated as a numeric column."""
    return name in _NUMFIELDS or name in _TIMEFIELDS or name in _ELEMENTS


def _parse_excel(raw: pd.DataFrame) -> pd.DataFrame:
    """Parse a transposed geochemical template DataFrame.

    *raw* is the DataFrame as returned by :func:`pd.read_excel` with
    ``header=None``.  Rows are variables, columns are samples.

    Returns a tidy DataFrame (samples × variables).
    """
    # Find the 'reference' row (start of the data table)
    start = None
    for i, val in enumerate(raw.iloc[:, 0]):
        if str(val).strip().lower() == 'reference':
            start = i
            break
    if start is None:
        raise ValueError("Could not find 'reference' row in Excel template.")

    block = raw.iloc[start:].reset_index(drop=True)

    n_samples = block.shape[1] - 1   # first col = var names

    cols: dict[str, list] = {}

    for _, row in block.iterrows():
        raw_name = str(row.iloc[0]).strip()
        if not raw_name or raw_name.lower() == 'nan':
            continue

        name = _clean_name(raw_name)
        name = _resolve_name(name)

        values = list(row.iloc[1:])   # one value per sample

        # Determine if this is a numeric field (after unit stripping)
        name_no_unit = re.sub(r'_(ppb|ppt|per|wtper|ppm)$', '', name)
        is_num = _is_numeric_field(name_no_unit)

        if is_num:
            # Apply BDL handling and '<' / '>' substitution
            cleaned = []
            for v in values:
                sv = str(v).strip().lower()
                if sv in _BDL_TOKENS:
                    cleaned.append(0.0)
                elif sv == '-' or sv == 'nan':
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
                        cleaned.append(np.nan)
            series = pd.Series(cleaned, dtype=float)
            # Apply unit conversion
            name, series = _apply_units(name, series)
            # Add _ppm suffix for elemental columns
            bare = re.sub(r'_(ppb|ppt|per|wtper|ppm)$', '', name)
            if bare in _ELEMENTS:
                name = bare + '_ppm'
            values = list(series)
        else:
            # Text field — keep as strings, drop pure NaN values
            values = [
                '' if str(v).strip().lower() in ('nan', '') else str(v).strip()
                for v in values
            ]

        if name in cols:
            # Duplicate: merge by averaging numeric, keep first for text
            existing = cols[name]
            if is_num:
                merged = []
                for a, b in zip(existing, values):
                    if np.isnan(float(a)) if isinstance(a, (int, float)) else True:
                        merged.append(b)
                    elif np.isnan(float(b)) if isinstance(b, (int, float)) else False:
                        merged.append(a)
                    else:
                        try:
                            merged.append((float(a) + float(b)) / 2)
                        except (TypeError, ValueError):
                            merged.append(a)
                cols[name] = merged
            # For text fields keep the first occurrence
        else:
            cols[name] = values

    df = pd.DataFrame(cols)

    # Add required columns if missing
    for req in _REQUIRED:
        if req not in df.columns:
            df[req] = ''

    return df


def pgcsv(input_path, output_dir=None):
    """Convert geochemical Excel template file(s) to database-ready CSV.

    Parameters
    ----------
    input_path : str or path-like
        Path to a single ``.xls`` / ``.xlsx`` file, or a directory
        containing such files.  Files with other extensions are skipped.
    output_dir : str or path-like, optional
        Directory to write CSV files.  Defaults to the same directory as
        each input file.

    Returns
    -------
    list of pathlib.Path
        Paths of the CSV files written.
    """
    input_path = pathlib.Path(input_path)

    if input_path.is_dir():
        files = sorted(
            f for f in input_path.iterdir()
            if f.suffix.lower() in ('.xls', '.xlsx')
        )
    elif input_path.is_file():
        files = [input_path]
    else:
        raise FileNotFoundError(f"Not found: {input_path}")

    written = []
    for f in files:
        try:
            raw = pd.read_excel(f, header=None, dtype=object)
            df = _parse_excel(raw)

            out_dir = pathlib.Path(output_dir) if output_dir else f.parent
            out_dir.mkdir(parents=True, exist_ok=True)
            out_path = out_dir / (f.stem + '.csv')
            df.to_csv(out_path, index=False, encoding='utf-8')
            written.append(out_path)
            print(f"Wrote {out_path}  ({len(df)} samples, {len(df.columns)} columns)")
        except Exception as exc:
            warnings.warn(f"Skipping {f.name}: {exc}")

    return written
