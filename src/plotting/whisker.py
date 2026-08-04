"""
Binned box-and-whisker plot utility.

Port of matlab/plotting/whisker.m.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


_DEFAULT_QUANTILES = [0.05, 0.25, 0.50, 0.75, 0.95]


def _gausscensor_quantiles(x_raw, quantiles):
    """Compute log₁₀ quantiles of left-censored data via gausscensor.

    Uses the Michael & Schucany (1986) empirical CDF method, weakly assuming
    log-normality.  Returns ``(parametric_q, empirical_q)`` both in log₁₀ space,
    or ``(None, None)`` if fewer than 8 detected values are present.

    Parameters
    ----------
    x_raw : array-like
        Values with BDL convention: positive = detected, negative = BDL with
        known detection limit (|value|), zero = BDL with unknown detection limit.
    quantiles : sequence of float
        Probability levels (ascending, matching whisker's 5-element array).
    """
    from src.geochem.processing.censored import gausscensor

    x = np.asarray(x_raw, dtype=float)
    x = x[np.isfinite(x)]
    if (x > 0).sum() < 8:
        return None, None

    try:
        _, qvals = gausscensor(x, scale='log10', q=list(quantiles))
        # qvals shape: (n_q, 2) — col 0 = parametric, col 1 = empirical (interpolated)
        empirical = qvals[:, 1]
        parametric = qvals[:, 0]
        if np.all(np.isfinite(empirical)):
            return parametric, empirical
        if np.all(np.isfinite(parametric)):
            return parametric, parametric
    except Exception:
        pass
    return None, None


def whisker(
    x_bins,
    y_bins,
    quantiles=None,
    scale='linear',
    color=None,
    plot=True,
    ax=None,
    y_raw=None,
):
    """Compute quantiles and draw 2-D box-and-whisker plots for binned data.

    Port of MATLAB ``whisker`` (cell-array form only).  For each bin,
    both the x and y distributions are summarised by five quantiles; a
    2-D box (IQR × IQR) with crossing whiskers (5th–95th %) is drawn
    at each bin's bivariate median.

    Parameters
    ----------
    x_bins : list of array-like
        One array of x-values per bin.
    y_bins : list of array-like
        One array of y-values per bin (same length as *x_bins*).
    quantiles : sequence of five floats, optional
        Ascending quantile levels.  Default is
        ``[0.05, 0.25, 0.50, 0.75, 0.95]``.
    scale : {'linear', 'log'}
        Apply log₁₀ to **y** values before computing quantiles.
    color : colour spec, optional
        Box / whisker colour.  Default: ``(0, 0.447, 0.741)``.
    plot : bool
        If False, compute quantiles only (no drawing).
    ax : matplotlib Axes, optional
    y_raw : list of array-like, optional
        Per-bin raw Y values with the BDL convention used throughout this
        package (positive = detected, negative = BDL with known detection
        limit, zero = BDL with unknown detection limit).  When provided
        **and** ``scale='log'``, quantiles are estimated via
        ``gausscensor`` (Michael & Schucany 1986) rather than
        ``numpy.quantile``.  Falls back to the standard detected-only
        quantile if fewer than 8 detected values are in the bin.

    Returns
    -------
    X : ndarray, shape (n_bins, 5)
        Quantiles of x in each bin.
    Y : ndarray, shape (n_bins, 5)
        Quantiles of y (or log₁₀(y)) in each bin.
    """
    if color is None:
        color = (0.0, 0.447, 0.741)
    if quantiles is None:
        quantiles = _DEFAULT_QUANTILES
    quantiles = np.asarray(quantiles, dtype=float)
    if quantiles.shape != (5,):
        raise ValueError("quantiles must have exactly 5 elements")

    use_censored = (y_raw is not None) and (scale == 'log')

    n_bins = len(x_bins)
    X = np.full((n_bins, 5), np.nan)
    Y = np.full((n_bins, 5), np.nan)

    for i, (xb, yb) in enumerate(zip(x_bins, y_bins)):
        xb = np.asarray(xb, dtype=float)
        yb = np.asarray(yb, dtype=float)

        if scale == 'log':
            yb = np.where(yb > 0, np.log10(yb), np.nan)

        mask = np.isfinite(xb) & np.isfinite(yb)
        if mask.sum() < 2:
            continue

        X[i] = np.quantile(xb[mask], quantiles)

        # Y quantiles: censored path or standard np.quantile
        y_set = False
        if use_censored:
            _, emp_q = _gausscensor_quantiles(y_raw[i], list(quantiles))
            if emp_q is not None:
                Y[i] = emp_q
                y_set = True

        if not y_set:
            Y[i] = np.quantile(yb[mask], quantiles)

    if not plot:
        return X, Y

    if ax is None:
        ax = plt.gca()

    light = np.clip(np.array(color) + 0.25, 0, 1)

    for i in range(n_bins):
        if np.any(np.isnan(X[i])) or np.any(np.isnan(Y[i])):
            continue

        # horizontal whisker at y-median spanning x[5%] – x[95%]
        ax.plot([X[i, 0], X[i, 4]], [Y[i, 2], Y[i, 2]], '-', color=color,
                linewidth=0.75)
        # vertical whisker at x-median spanning y[5%] – y[95%]
        ax.plot([X[i, 2], X[i, 2]], [Y[i, 0], Y[i, 4]], '-', color=color,
                linewidth=0.75)
        # IQR box (x[25%:75%] × y[25%:75%])
        box = Rectangle(
            (X[i, 1], Y[i, 1]),
            X[i, 3] - X[i, 1],
            Y[i, 3] - Y[i, 1],
            facecolor=light, edgecolor=color, linewidth=0.75, alpha=0.6,
        )
        ax.add_patch(box)
        # median point
        ax.plot(X[i, 2], Y[i, 2], 'o',
                markerfacecolor='white', markeredgecolor=color,
                markersize=4, linewidth=0.75)

    return X, Y
