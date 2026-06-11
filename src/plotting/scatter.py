"""
Scatter plot helpers.

:func:`scatterbar`
    Scatter plot with optional symmetric or asymmetric error bars in x,
    y, or both.  Equivalent to MATLAB ``scatterbar.m``.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def scatterbar(
    x, y,
    xerr=None, yerr=None,
    xerr_upper=None, yerr_upper=None,
    err_type='sigma',
    color=None,
    size=None,
    marker='o',
    bar_color=(0.8, 0.8, 0.8),
    bar_lw=0.8,
    ax=None,
    **scatter_kwargs,
):
    """Scatter plot with optional x and/or y error bars.

    Error bars are drawn first (in grey) and the scatter markers on top,
    matching the MATLAB ``scatterbar.m`` layering.

    Parameters
    ----------
    x, y : array-like
        Data coordinates.
    xerr, yerr : array-like or None
        Symmetric error (one value per point) or lower bound when
        *xerr_upper* / *yerr_upper* is also provided.
    xerr_upper, yerr_upper : array-like or None
        Upper bound for asymmetric errors.  Ignored when *xerr* / *yerr*
        is None.
    err_type : {'sigma', 'value'}
        ``'sigma'`` (default): errors are half-widths, so the bar extends
        from ``v - err_lower`` to ``v + err_upper``.
        ``'value'``: errors are the actual endpoint coordinates.
    color : array-like, optional
        Marker colour(s).  Accepts anything matplotlib scatter accepts.
    size : scalar or array-like, optional
        Marker size(s) in points² (default: matplotlib default).
    marker : str
        Marker style (default ``'o'``).
    bar_color : color-like
        Error bar line colour (default light grey).
    bar_lw : float
        Error bar line width (default 0.8).
    ax : matplotlib.axes.Axes, optional
        Target axes.  Uses current axes if None.
    **scatter_kwargs
        Additional keyword arguments forwarded to ``ax.scatter()``.

    Returns
    -------
    scatter_handle : matplotlib.collections.PathCollection
    """
    if ax is None:
        ax = plt.gca()

    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()

    bar_kw = dict(color=bar_color, linewidth=bar_lw, zorder=1)

    # ---- x error bars ----
    if xerr is not None:
        xl, xu = _resolve_errors(x, xerr, xerr_upper, err_type)
        for i in range(len(x)):
            if np.isfinite(xl[i]) and np.isfinite(xu[i]):
                ax.plot([xl[i], xu[i]], [y[i], y[i]], **bar_kw)

    # ---- y error bars ----
    if yerr is not None:
        yl, yu = _resolve_errors(y, yerr, yerr_upper, err_type)
        for i in range(len(y)):
            if np.isfinite(yl[i]) and np.isfinite(yu[i]):
                ax.plot([x[i], x[i]], [yl[i], yu[i]], **bar_kw)

    # ---- scatter markers ----
    sc_kw = dict(zorder=2, marker=marker)
    if size is not None:
        sc_kw['s'] = size
    if color is not None:
        sc_kw['c'] = color
    sc_kw.update(scatter_kwargs)

    handle = ax.scatter(x, y, **sc_kw)
    return handle


def _resolve_errors(v, err_lower, err_upper, err_type):
    """Return (lower_endpoint, upper_endpoint) arrays."""
    err_lower = np.asarray(err_lower, dtype=float).ravel()
    if err_upper is None:
        err_upper = err_lower.copy()
    else:
        err_upper = np.asarray(err_upper, dtype=float).ravel()

    if err_type == 'sigma':
        return v - err_lower, v + err_upper
    elif err_type == 'value':
        return err_lower, err_upper
    else:
        raise ValueError(f"Unknown err_type '{err_type}'.  Use 'sigma' or 'value'.")
