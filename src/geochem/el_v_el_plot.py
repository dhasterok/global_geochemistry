"""
Element-vs-element (and ratio-vs-element) plots.

Main entry point:

    ax = elvelplot(data, 'sio2', 'mgo')
    ax = elvelplot(data, ['la', 'yb'], 'sio2', plotstyle='whisker',
                   xscale='log', binfield='rock_group')

``el1`` / ``el2`` can be a column name (str) or a two-element list/tuple to
form a ratio: ``['la', 'yb']`` → La/Yb plotted on the axis.

Four plot styles
----------------
``'scatter'``
    Filled semi-transparent scatter of all detected (> 0) values.
``'hist2d'``
    2-D density heatmap.  Bins are computed by the Freedman–Diaconis rule
    on the data in the axis-scale space.
``'average'``
    Median point with IQR crosshairs (grey lines), one point per bin.
    Quantiles estimated via :func:`~.processing.censored.gausscensor` to
    account for below-detection-limit values.
``'whisker'``
    2.5–97.5 % extent lines (grey), filled 25–75 % box, median white dot.
    Same censored-quantile estimation as ``'average'``.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.patches import Rectangle

from global_geochemistry.geochem.processing.censored import gausscensor


# --------------------------------------------------------------------------
# Label helper (mirrors MATLAB getlabel.m)
# --------------------------------------------------------------------------

_OXIDE_COLS = {'sio2', 'tio2', 'al2o3', 'feo_tot', 'mgo', 'cao', 'mno',
               'na2o', 'k2o', 'p2o5', 'loi'}

_LABEL_MAP = {
    'heat_production':  'Heat Production (μW m⁻³)',
    'density_model':    'Density (kg m⁻³)',
    'crustal_thickness':'Crustal Thickness (km)',
    'p_velocity':       'P Velocity (km/s)',
    's_velocity':       'S Velocity (km/s)',
    'mg_number':        'Mg Number',
}

def _make_label(el) -> str:
    if isinstance(el, (list, tuple)):
        return f'{el[0]}/{el[1]}'
    col = str(el).lower()
    if col in _LABEL_MAP:
        return _LABEL_MAP[col]
    if col in _OXIDE_COLS:
        return f'{el} (wt.%)'
    return f'{el} (ppm)'


# --------------------------------------------------------------------------
# Data extraction helpers
# --------------------------------------------------------------------------

def _extract(data: pd.DataFrame, el, transform_log: bool) -> np.ndarray:
    """Return 1-D array of detected values (>0), optionally log10-transformed."""
    if isinstance(el, (list, tuple)):
        a = data[el[0]].to_numpy(float)
        b = data[el[1]].to_numpy(float)
        mask = (a > 0) & (b > 0)
        v = np.where(mask, a / b, np.nan)
    else:
        raw = data[el].to_numpy(float)
        v = np.where(raw > 0, raw, np.nan)

    if transform_log:
        with np.errstate(divide='ignore', invalid='ignore'):
            v = np.log10(v)
    return v


def _censored_quantiles(
    data: pd.DataFrame,
    el,
    scale: str,
    q=(0.025, 0.25, 0.5, 0.75, 0.975),
) -> np.ndarray | None:
    """Return shape (5,) quantile array, or None if not enough data.

    Values in the same space as *scale*: log10 for 'log', linear for
    'linear'.  Caller is responsible for back-transforming if needed.
    """
    if isinstance(el, (list, tuple)):
        a = data[el[0]].to_numpy(float)
        b = data[el[1]].to_numpy(float)
        mask = (a > 0) & (b > 0)
        v = (a / b)[mask]
    else:
        raw = data[el].to_numpy(float)
        # keep all non-NaN values; gausscensor handles BDL (<=0) internally
        v = raw[~np.isnan(raw)]

    censor_scale = 'log10' if scale == 'log' else 'linear'
    result = gausscensor(v, scale=censor_scale, q=list(q))
    if not isinstance(result, tuple):
        return None

    _, qvals = result
    # qvals is (N_q, 2); column 1 = empirical, column 0 = Gaussian model
    quantiles = qvals[:, 1]

    # back-transform from log10 to linear so caller can use ax.set_xscale('log')
    if scale == 'log':
        quantiles = 10.0 ** quantiles
    return quantiles


# --------------------------------------------------------------------------
# Bin-edge helper (Freedman–Diaconis rule, mirrors MATLAB binedges.m)
# --------------------------------------------------------------------------

def _bin_edges_fd(v: np.ndarray) -> np.ndarray:
    """Freedman–Diaconis bin edges clipped to [p2.5, p97.5]."""
    v = v[np.isfinite(v)]
    if len(v) < 4:
        return np.linspace(float(np.nanmin(v)), float(np.nanmax(v)), 11)
    qs = np.quantile(v, [0.025, 0.25, 0.75, 0.975])
    vmin, vmax = np.floor(qs[0]), np.ceil(qs[3])
    iqr = qs[2] - qs[1]
    dv = 2.0 * iqr * len(v) ** (-1.0 / 3.0)
    if dv <= 0 or not np.isfinite(dv):
        return np.linspace(vmin, vmax, 21)
    nv = max(int(np.ceil((vmax - vmin) / dv)), 2)
    return np.linspace(vmin, vmax, nv)


# --------------------------------------------------------------------------
# Main function
# --------------------------------------------------------------------------

_QUANTS = (0.025, 0.25, 0.5, 0.75, 0.975)


def elvelplot(
    data: pd.DataFrame,
    el1,
    el2,
    plotstyle: str = 'scatter',
    xscale: str = 'linear',
    yscale: str = 'linear',
    binfield: str | None = None,
    bin_edges=None,
    color=None,
    min_n: int = 10,
    ax=None,
    **scatter_kw,
):
    """Plot element (or ratio) vs. element (or ratio).

    Parameters
    ----------
    data : pandas.DataFrame
    el1, el2 : str or list of two str
        Column name(s).  Two-element list → ratio ``el1[0]/el1[1]``.
    plotstyle : {'scatter', 'hist2d', 'average', 'whisker'}
    xscale, yscale : {'linear', 'log'}
    binfield : str, optional
        Column used to split data into bins for ``'average'`` / ``'whisker'``
        modes.  Numeric columns use *bin_edges* for binning; string/category
        columns bin by unique value.
    bin_edges : array-like, optional
        Bin edges for numeric *binfield*.  Defaults to midpoints of a
        Freedman–Diaconis grid.
    color : scalar, array-like, or matplotlib color, optional
        For ``'scatter'``: per-point values (colormap) or single colour.
        For ``'average'``/``'whisker'``: per-bin scalar values (colormap) or
        single colour.
    min_n : int
        Minimum samples per bin for stat modes (default 10).
    ax : matplotlib.axes.Axes, optional
    **scatter_kw
        Extra keyword arguments forwarded to :func:`matplotlib.axes.Axes.scatter`.

    Returns
    -------
    ax : matplotlib.axes.Axes
    stats : dict or None
        For ``'average'`` and ``'whisker'`` modes: ``{'xq': (5, N), 'yq': (5, N), 'n': (N,)}``.
        ``None`` for scatter and hist2d.

    Examples
    --------
    >>> fig, ax = plt.subplots()
    >>> elvelplot(data, 'sio2', 'mgo', ax=ax)

    >>> ax, stats = elvelplot(data, 'th', 'u',
    ...                       plotstyle='whisker', xscale='log', yscale='log',
    ...                       binfield='rock_group')
    """
    if ax is None:
        _, ax = plt.subplots()

    plotstyle = plotstyle.lower()

    if plotstyle in ('scatter', 'hist2d'):
        _plot_all(ax, data, el1, el2, xscale, yscale, plotstyle, color, scatter_kw)
        ax.set_xlabel(_make_label(el1))
        ax.set_ylabel(_make_label(el2))
        ax.set_box_aspect(None)
        if plotstyle == 'scatter':
            ax.set_xscale(xscale)
            ax.set_yscale(yscale)
        else:
            _set_log_ticks(ax, xscale, yscale)
        return ax, None

    elif plotstyle in ('average', 'whisker'):
        xq, yq, n = _plot_avg(ax, data, el1, el2, xscale, yscale,
                               plotstyle, binfield, bin_edges, color, min_n)
        ax.set_xlabel(_make_label(el1))
        ax.set_ylabel(_make_label(el2))
        ax.set_xscale(xscale)
        ax.set_yscale(yscale)
        ax.set_box_aspect(None)
        stats = {'xq': xq, 'yq': yq, 'n': n}
        return ax, stats

    else:
        raise ValueError(f"Unknown plotstyle '{plotstyle}'. "
                         "Choose 'scatter', 'hist2d', 'average', or 'whisker'.")


# --------------------------------------------------------------------------
# Scatter / hist2d
# --------------------------------------------------------------------------

def _plot_all(ax, data, el1, el2, xscale, yscale, plotstyle, color, scatter_kw):
    log_x = xscale == 'log'
    log_y = yscale == 'log'

    if plotstyle == 'scatter':
        x_raw = _extract(data, el1, transform_log=False)
        y_raw = _extract(data, el2, transform_log=False)
        mask = np.isfinite(x_raw) & np.isfinite(y_raw)
        x, y = x_raw[mask], y_raw[mask]

        kw = dict(s=18, linewidths=0, alpha=0.5)
        kw.update(scatter_kw)

        show_cbar = False
        if color is None:
            kw.setdefault('c', [[0.5, 0.5, 0.5]] * len(x))
        else:
            c = np.asarray(color)
            if c.ndim == 1 and len(c) == len(x_raw):
                # per-point scalar values → colormap
                kw['c'] = c[mask]
                show_cbar = True
            elif c.ndim == 2 and c.shape[0] == len(x_raw):
                # per-point RGB
                kw['c'] = c[mask]
            else:
                # single colour spec
                kw['c'] = [color] * len(x)

        sc = ax.scatter(x, y, **kw)
        if show_cbar:
            plt.colorbar(sc, ax=ax)

    else:  # hist2d
        x = _extract(data, el1, transform_log=log_x)
        y = _extract(data, el2, transform_log=log_y)
        mask = np.isfinite(x) & np.isfinite(y)
        x, y = x[mask], y[mask]
        if len(x) < 4:
            return

        ex = _bin_edges_fd(x)
        ey = _bin_edges_fd(y)
        counts, xedges, yedges = np.histogram2d(x, y, bins=[ex, ey])

        img = counts.T  # (n_ybins, n_xbins) → imshow rows=y, cols=x
        with np.errstate(divide='ignore'):
            img_log = np.where(img > 0, np.log10(img), np.nan)

        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
        ax.imshow(img_log, origin='lower', extent=extent, aspect='auto',
                  interpolation='nearest', cmap='plasma')
        ax.set_xlim(xedges[0], xedges[-1])
        ax.set_ylim(yedges[0], yedges[-1])


def _set_log_ticks(ax, xscale, yscale):
    """Replace log10-value tick labels with 10^n for hist2d on log scale."""
    if xscale == 'log':
        xl = ax.get_xlim()
        ticks = np.arange(np.ceil(xl[0]), np.floor(xl[1]) + 1, dtype=int)
        ax.set_xticks(ticks)
        ax.xaxis.set_major_formatter(
            mticker.FuncFormatter(lambda v, _: f'$10^{{{int(v)}}}$')
        )
    if yscale == 'log':
        yl = ax.get_ylim()
        ticks = np.arange(np.ceil(yl[0]), np.floor(yl[1]) + 1, dtype=int)
        ax.set_yticks(ticks)
        ax.yaxis.set_major_formatter(
            mticker.FuncFormatter(lambda v, _: f'$10^{{{int(v)}}}$')
        )


# --------------------------------------------------------------------------
# Average / whisker
# --------------------------------------------------------------------------

def _compute_bins(data, binfield, bin_edges):
    """Return list of (mask, bin_centre_or_label) tuples."""
    col = data[binfield]

    # Categorical
    if col.dtype == object or str(col.dtype) == 'category':
        bins = []
        for val in col.dropna().unique():
            mask = col == val
            if mask.sum() > 0:
                bins.append((mask.to_numpy(bool), val))
        return bins

    # Numeric
    vals = col.to_numpy(float)
    if bin_edges is None:
        finite = vals[np.isfinite(vals)]
        if len(finite) < 4:
            return []
        bin_edges = _bin_edges_fd(finite)

    edges = np.asarray(bin_edges, dtype=float)
    bins = []
    for i in range(len(edges) - 1):
        mask = (vals >= edges[i]) & (vals < edges[i + 1])
        centre = 0.5 * (edges[i] + edges[i + 1])
        bins.append((mask, centre))
    return bins


def _plot_avg(ax, data, el1, el2, xscale, yscale, plotstyle,
              binfield, bin_edges, color, min_n):

    if binfield is None:
        # Single set of quantiles for the whole dataset
        xq_full = _censored_quantiles(data, el1, xscale)
        yq_full = _censored_quantiles(data, el2, yscale)
        if xq_full is None or yq_full is None:
            return np.empty((5, 0)), np.empty((5, 0)), np.array([])
        xq = xq_full[:, np.newaxis]   # (5, 1)
        yq = yq_full[:, np.newaxis]
        n  = np.array([len(data)])
        c_vals = [color if color is not None else [0.2, 0.5, 0.8]]
        _draw_stat_glyphs(ax, xq, yq, c_vals, plotstyle)
        return xq, yq, n

    # Binned
    bins = _compute_bins(data, binfield, bin_edges)
    xqs, yqs, ns, centres = [], [], [], []

    for mask, centre in bins:
        subset = data[mask]
        if mask.sum() < min_n:
            continue
        xq_i = _censored_quantiles(subset, el1, xscale)
        yq_i = _censored_quantiles(subset, el2, yscale)
        if xq_i is None or yq_i is None:
            continue
        xqs.append(xq_i)
        yqs.append(yq_i)
        ns.append(mask.sum())
        centres.append(centre)

    if not xqs:
        return np.empty((5, 0)), np.empty((5, 0)), np.array([])

    xq = np.column_stack(xqs)   # (5, N_bins)
    yq = np.column_stack(yqs)
    n  = np.array(ns)

    # Colour mapping
    if color is not None:
        c_vals = [color] * len(centres)
    else:
        numeric_centres = [c for c in centres if isinstance(c, (int, float, np.floating))]
        if numeric_centres and len(numeric_centres) == len(centres):
            cmap = plt.colormaps['viridis']
            norm = plt.Normalize(min(centres), max(centres))
            c_vals = [cmap(norm(c)) for c in centres]
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            plt.colorbar(sm, ax=ax)
        else:
            cmap = plt.colormaps['tab10']
            c_vals = [cmap(i % 10) for i in range(len(centres))]

    _draw_stat_glyphs(ax, xq, yq, c_vals, plotstyle)
    return xq, yq, n


def _draw_stat_glyphs(ax, xq, yq, c_vals, plotstyle):
    """Draw either 'average' crosshairs or 'whisker' boxes onto *ax*.

    xq, yq : (5, N_bins) — rows are [p2.5, p25, p50, p75, p97.5]
    c_vals  : list of N_bins colours
    """
    grey = [0.7, 0.7, 0.7]

    for i in range(xq.shape[1]):
        x = xq[:, i]
        y = yq[:, i]
        c = c_vals[i] if i < len(c_vals) else grey

        if plotstyle == 'average':
            # IQR crosshairs in grey
            ax.plot([x[1], x[3]], [y[2], y[2]], '-', color=grey, linewidth=0.75, zorder=1)
            ax.plot([x[2], x[2]], [y[1], y[3]], '-', color=grey, linewidth=0.75, zorder=1)
            ax.scatter(x[2], y[2], c=[c], s=40, zorder=3, linewidths=0)

        else:  # whisker
            # 2.5–97.5 % extent lines
            ax.plot([x[0], x[4]], [y[2], y[2]], '-', color=grey, linewidth=0.75, zorder=1)
            ax.plot([x[2], x[2]], [y[0], y[4]], '-', color=grey, linewidth=0.75, zorder=1)
            # 25–75 % filled box — use fill (not Rectangle) for log-axis compatibility
            bx = [x[1], x[3], x[3], x[1], x[1]]
            by = [y[1], y[1], y[3], y[3], y[1]]
            ax.fill(bx, by, color=[0.2, 0.5, 0.8], zorder=2, linewidth=0)
            # Median white dot
            ax.scatter(x[2], y[2], c=[[1, 1, 1]], s=25, zorder=4,
                       edgecolors=grey, linewidths=0.5)
