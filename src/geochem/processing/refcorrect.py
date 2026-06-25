"""
Reference-variable corrections for geochemical concentrations.

Two approaches, both generalised to any reference field (not just SiO₂):

``poly_correct(data, el, ref, target, ...)``
    Fits a weighted polynomial to binned medians of the response vs. the
    reference variable and shifts each sample by the fitted trend difference:

        y_corrected = y × 10^( f(target) − f(x) )     [log_y=True]
        y_corrected = y  +  f(target) − f(x)           [log_y=False]

    Fast and interpretable; the polynomial coefficients fully describe the
    correction.  Equivalent to (and generalises) MATLAB ``sio2_correction.m``.

``quantile_correct(data, elements, ref, target, ...)``
    Rank-preserving correction.  For each element, fits smooth polynomial
    models for the log-normal location μ(x) and scale σ(x) parameters across
    bins of the reference variable (using :func:`~.censored.gausscensor` to
    handle below-detection-limit values).  Each sample is then mapped through
    its z-score to the target reference value:

        z = ( log₁₀(y) − μ(x) ) / σ(x)
        y_corrected = 10^( μ(x_target) + σ(x_target) × z )

    Preserves each sample's rank within the local population.  Equivalent to
    (and generalises) MATLAB ``tracenorm.m`` / ``adjbyquantile2.m``.

Both functions:

* Accept any reference column (SiO₂, MgO, age, …)
* Support optional log₁₀ transforms on the response (``log_y``) and/or the
  reference variable (``log_x``)
* Support arbitrary polynomial degree for the trend fitting
* Weight bin fits by √N to down-weight sparse bins
* Return corrected values without producing plots (plotting is separate)
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from global_geochemistry.geochem.processing.censored import gausscensor


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _to_log(v: np.ndarray) -> np.ndarray:
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(v > 0, np.log10(v), np.nan)


def _from_log(v: np.ndarray) -> np.ndarray:
    return 10.0 ** v


def _uniform_bins(x_valid: np.ndarray, n_bins: int):
    """Return (edges, midpoints) for n_bins uniform bins over x_valid range."""
    edges = np.linspace(x_valid.min(), x_valid.max(), n_bins + 1)
    mids  = 0.5 * (edges[:-1] + edges[1:])
    return edges, mids


def _wpolyfit(x: np.ndarray, y: np.ndarray, degree: int,
              w: np.ndarray | None = None) -> np.ndarray:
    """Weighted polynomial fit; returns coefficients in np.polyval convention."""
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() <= degree:
        return np.full(degree + 1, np.nan)
    xm, ym = x[mask], y[mask]
    wm = w[mask] if w is not None else np.ones_like(xm)
    return np.polyfit(xm, ym, degree, w=wm)


# ---------------------------------------------------------------------------
# poly_correct
# ---------------------------------------------------------------------------

def poly_correct(
    data: pd.DataFrame,
    el: str,
    ref: str,
    target: float,
    degree: int = 2,
    log_y: bool = True,
    log_x: bool = False,
    n_bins: int = 20,
    min_n: int = 5,
) -> tuple[np.ndarray, np.ndarray]:
    """Correct concentrations by removing a polynomial trend with a reference variable.

    Bins the data by *ref*, computes the median of *el* (or log₁₀(*el*) when
    *log_y* is ``True``) in each bin, fits a weighted polynomial of degree
    *degree* to those bin medians, then shifts every sample to *target*:

        y_corrected = y × 10^( f(target) − f(x) )   [log_y=True]
        y_corrected = y  +  f(target) − f(x)          [log_y=False]

    Parameters
    ----------
    data : pandas.DataFrame
    el : str
        Column to correct (response variable).
    ref : str
        Reference column (e.g. ``'sio2'``, ``'mgo'``, ``'avg_age'``).
    target : float
        Value of *ref* to normalise to (in the same units / scale as *ref*,
        *before* any log₁₀ transform).
    degree : int
        Polynomial degree (1 = linear, 2 = quadratic, …).  Default 2.
    log_y : bool
        Fit and correct in log₁₀(*el*) space (default ``True``).
        Samples with *el* ≤ 0 are treated as NaN.
    log_x : bool
        Apply log₁₀ to *ref* before binning and fitting (default ``False``).
        Useful when the reference variable has a log-normal distribution.
        *target* is always supplied in original (linear) units.
    n_bins : int
        Number of uniform bins along *ref* for computing bin medians (default 20).
    min_n : int
        Minimum samples per bin for inclusion in the polynomial fit (default 5).

    Returns
    -------
    y_corrected : numpy.ndarray
        Corrected values in the original (linear) units of *el*.  NaN where
        the input was NaN or ≤ 0 (when *log_y* is ``True``).
    poly_coeffs : numpy.ndarray
        Polynomial coefficients in :func:`numpy.polyval` convention
        (highest degree first).

    Examples
    --------
    >>> hp_corr, coeffs = poly_correct(data, 'heat_production', 'sio2', target=60.0)
    >>> th_corr, _      = poly_correct(data, 'th_ppm', 'mgo', target=5.0, degree=1)
    >>> u_corr, _       = poly_correct(data, 'u_ppm', 'avg_age', target=500.0,
    ...                                log_x=True)
    """
    y_raw = data[el].to_numpy(float)
    x_raw = data[ref].to_numpy(float)

    # Working copies in transformed space
    y_work = _to_log(y_raw) if log_y else y_raw.copy()
    x_work = _to_log(x_raw) if log_x else x_raw.copy()
    x_tgt  = np.log10(target) if log_x else float(target)

    valid  = np.isfinite(y_work) & np.isfinite(x_work)

    if valid.sum() <= degree + 1:
        return y_raw.copy(), np.full(degree + 1, np.nan)

    edges, mids = _uniform_bins(x_work[valid], n_bins)

    # Per-bin medians and sample counts
    y_bin = np.full(n_bins, np.nan)
    n_bin = np.zeros(n_bins, dtype=float)
    for i in range(n_bins):
        mask = (edges[i] <= x_work) & (x_work < edges[i + 1]) & valid
        n_bin[i] = mask.sum()
        if n_bin[i] >= min_n:
            y_bin[i] = np.nanmedian(y_work[mask])

    ok = np.isfinite(y_bin) & (n_bin >= min_n)
    coeffs = _wpolyfit(mids[ok], y_bin[ok], degree, w=np.sqrt(n_bin[ok]))

    # Shift each sample
    f_target = np.polyval(coeffs, x_tgt)
    f_x      = np.polyval(coeffs, x_work)
    delta    = f_target - f_x  # in y_work space

    y_corr = y_work.copy()
    y_corr[valid] += delta[valid]

    y_corrected = _from_log(y_corr) if log_y else y_corr
    # Restore NaN / negative inputs
    y_corrected[~valid] = y_raw[~valid]

    return y_corrected, coeffs


# ---------------------------------------------------------------------------
# quantile_correct
# ---------------------------------------------------------------------------

def quantile_correct(
    data: pd.DataFrame,
    elements: list[str] | str,
    ref: str,
    target: float,
    mu_degree: int = 2,
    sigma_degree: int = 1,
    log_y: bool = True,
    log_x: bool = False,
    n_bins: int = 20,
    min_n: int = 8,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Rank-preserving correction of concentrations to a reference value.

    For each element in *elements*, fits smooth polynomial models for the
    log-normal location μ(x) and scale σ(x) using
    :func:`~.censored.gausscensor` per bin of the reference variable.  Each
    sample is then mapped through its z-score to the target:

        z = ( log₁₀(y) − μ(x) ) / σ(x)
        y_corrected = 10^( μ(target) + σ(target) × z )

    This preserves the sample's rank within the local population and accounts
    for changes in both the typical value *and* spread across the reference
    variable.  Equivalent to (and generalises) MATLAB ``adjbyquantile2``.

    Parameters
    ----------
    data : pandas.DataFrame
    elements : str or list of str
        Column(s) to correct.
    ref : str
        Reference column (e.g. ``'sio2'``, ``'mgo'``).
    target : float
        Value of *ref* to normalise to (linear units).
    mu_degree : int
        Polynomial degree for the μ(x) model (default 2 = quadratic).
    sigma_degree : int
        Polynomial degree for the σ(x) model (default 1 = linear).
    log_y : bool
        Work in log₁₀(*el*) space (default ``True``).  When ``False``, the
        distribution is assumed to be Gaussian in linear space (rarely
        appropriate for trace elements).
    log_x : bool
        Apply log₁₀ to *ref* before fitting (default ``False``).
    n_bins : int
        Number of uniform bins (default 20).
    min_n : int
        Minimum samples per bin for gausscensor (default 8; fewer than ~8
        points makes the censored-data CDF estimate unreliable).

    Returns
    -------
    corrected : pandas.DataFrame
        Copy of *data* with the corrected element columns.  Non-element
        columns are preserved unchanged.  Samples outside the *ref* range
        or with NaN *ref* values are left as NaN.
    models : pandas.DataFrame
        One row per element with columns:
        ``element``, ``mu_coeffs`` (array), ``sigma_coeffs`` (array),
        ``mu_rms``, ``sigma_rms``, ``n_bins_used``.

    Examples
    --------
    >>> corr, models = quantile_correct(data,
    ...     ['th_ppm', 'u_ppm', 'la_ppm', 'yb_ppm'],
    ...     ref='sio2', target=60.0)
    >>> corr, _ = quantile_correct(data, 'heat_production', ref='mgo',
    ...                             target=5.0, mu_degree=3)
    """
    if isinstance(elements, str):
        elements = [elements]

    x_raw  = data[ref].to_numpy(float)
    x_work = _to_log(x_raw) if log_x else x_raw.copy()
    x_tgt  = np.log10(target) if log_x else float(target)

    corrected = data.copy()
    model_rows: list[dict] = []

    for el in elements:
        y_raw  = data[el].to_numpy(float)
        y_work = _to_log(y_raw) if log_y else y_raw.copy()

        valid  = np.isfinite(y_work) & np.isfinite(x_work)

        if valid.sum() <= max(mu_degree, sigma_degree) + 1:
            model_rows.append({
                'element': el,
                'mu_coeffs': np.full(mu_degree + 1, np.nan),
                'sigma_coeffs': np.full(sigma_degree + 1, np.nan),
                'mu_rms': np.nan, 'sigma_rms': np.nan, 'n_bins_used': 0,
            })
            continue

        edges, mids = _uniform_bins(x_work[valid], n_bins)

        # Step 1: estimate log-normal parameters per bin via gausscensor
        mu_bin    = np.full(n_bins, np.nan)
        sigma_bin = np.full(n_bins, np.nan)
        n_bin     = np.zeros(n_bins, dtype=float)

        censor_scale = 'log10' if log_y else 'linear'

        for i in range(n_bins):
            mask = (edges[i] <= x_work) & (x_work < edges[i + 1]) & np.isfinite(x_work)
            n_bin[i] = mask.sum()
            if n_bin[i] < min_n:
                continue
            y_slice = y_raw[mask]  # pass original values; gausscensor handles BDL
            result = gausscensor(y_slice, scale=censor_scale)
            if not isinstance(result, tuple):
                continue
            model_i, _ = result
            if not np.isfinite(model_i.get('mu', np.nan)):
                continue
            mu_bin[i]    = model_i['mu']
            sigma_bin[i] = model_i['sigma']

        ok = np.isfinite(mu_bin) & np.isfinite(sigma_bin) & (n_bin >= min_n)
        n_used = ok.sum()

        w = np.sqrt(n_bin[ok])

        # Step 2: fit smooth polynomials to μ(x) and σ(x)
        mu_coeffs    = _wpolyfit(mids[ok], mu_bin[ok],    mu_degree,    w=w)
        sigma_coeffs = _wpolyfit(mids[ok], sigma_bin[ok], sigma_degree, w=w)

        # RMS of polynomial fits (in y_work / log10 units)
        mu_rms = sigma_rms = np.nan
        if n_used > mu_degree + 1:
            mu_pred = np.polyval(mu_coeffs, mids[ok])
            mu_rms  = float(np.sqrt(np.mean((mu_bin[ok] - mu_pred) ** 2)))
        if n_used > sigma_degree + 1:
            sg_pred   = np.polyval(sigma_coeffs, mids[ok])
            sigma_rms = float(np.sqrt(np.mean((sigma_bin[ok] - sg_pred) ** 2)))

        # Step 3: map each sample through its z-score to x_target
        mu_x     = np.polyval(mu_coeffs,    x_work)
        sigma_x  = np.polyval(sigma_coeffs, x_work)
        mu_tgt   = np.polyval(mu_coeffs,    x_tgt)
        sigma_tgt= np.polyval(sigma_coeffs, x_tgt)

        # z-score in transformed y space (log10 if log_y)
        with np.errstate(divide='ignore', invalid='ignore'):
            z = (y_work - mu_x) / sigma_x

        y_corr_work = mu_tgt + sigma_tgt * z

        if log_y:
            y_corrected = _from_log(y_corr_work)
        else:
            y_corrected = y_corr_work

        # Preserve BDL sign convention (negative → negative corrected)
        bdl = (y_raw < 0) & np.isfinite(y_raw)
        if bdl.any():
            # BDL: correct the absolute value, reapply sign
            y_bdl  = _to_log(-y_raw[bdl]) if log_y else -y_raw[bdl]
            z_bdl  = (y_bdl - mu_x[bdl]) / sigma_x[bdl]
            yc_bdl = mu_tgt + sigma_tgt * z_bdl
            y_corrected[bdl] = -(_from_log(yc_bdl) if log_y else yc_bdl)

        # NaN out samples outside the fitted range or with invalid ref
        outside = ~np.isfinite(x_work) | (x_work < edges[0]) | (x_work > edges[-1])
        y_corrected[outside | ~valid] = np.nan

        corrected[el] = y_corrected

        model_rows.append({
            'element': el,
            'mu_coeffs': mu_coeffs,
            'sigma_coeffs': sigma_coeffs,
            'mu_rms': mu_rms,
            'sigma_rms': sigma_rms,
            'n_bins_used': int(n_used),
        })

    models = pd.DataFrame(model_rows)
    return corrected, models
