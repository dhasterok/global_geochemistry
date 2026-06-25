"""
Radar (spider) and sequential radar plots.

Two functions:
- ``radar_prep``  —  aggregate a DataFrame into per-group stats for ``radarplot``
- ``radarplot``   —  draw the plot on a plain matplotlib Axes (no polar projection)

Two styles (``style`` parameter):
- ``'web'``     Classic spider/radar web.  Each group becomes a closed polygon.
                Quantile bands are drawn as shaded rings between inner/outer lines.
- ``'series'``  Sequential radar.  Each row of the input array is one lap around
                the spokes; consecutive laps are connected, with a colormap
                encoding the sequence so you can see how the pattern evolves.
                Direction (clockwise / counterclockwise) is selectable.

Examples
--------
Basic web from a DataFrame::

    import pandas as pd
    from src.plotting.radar import radar_prep, radarplot

    fields = ['SiO2', 'TiO2', 'Al2O3', 'FeOT', 'MgO', 'CaO', 'Na2O', 'K2O']
    prep = radar_prep(df, fields, group_field='rock_type',
                      groups=['basalt', 'andesite', 'rhyolite'],
                      quantiles=[0.25, 0.50, 0.75])
    fig, ax = plt.subplots(figsize=(7, 7))
    radarplot(**prep, ax=ax)

Set per-axis ranges and log scale::

    radarplot(**prep, ax=ax,
              axis_lim={'SiO2': (40, 80), 'TiO2': (0.1, 3)},
              log=['TiO2', 'K2O'],
              axis_labels={'TiO2': 'TiO₂ (log)'})

Sequential / time-series style::

    # data: 2-D array (n_obs, n_fields)
    fig, ax = plt.subplots(figsize=(7, 7))
    radarplot(vals=data, fields=fields, style='series',
              clockwise=True, cmap='plasma', ax=ax)
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import LineCollection
import matplotlib.colors as mcolors


# ---------------------------------------------------------------------------
# Data preparation
# ---------------------------------------------------------------------------

def radar_prep(
    data,
    fields: list[str],
    group_field: str | None = None,
    groups=None,
    quantiles: list[float] | None = None,
) -> dict:
    """Aggregate a DataFrame into per-group statistics for :func:`radarplot`.

    Parameters
    ----------
    data : pandas.DataFrame
        Source data table.
    fields : list of str
        Column names that form the spokes of the radar.
    group_field : str, optional
        Column used to partition rows into groups.  If *None*, the whole
        dataset is treated as a single group.
    groups : list, optional
        Values of *group_field* to include, in plotting order.  If *None*,
        all unique values are used (sorted).
    quantiles : list of float, optional
        Quantile levels, e.g. ``[0.25, 0.50, 0.75]``.  If *None*, the mean
        is computed.  Supported counts: 1, 2, 3, 5.

    Returns
    -------
    dict
        Keys ``'vals'``, ``'fields'``, ``'groups'``, ``'quantiles'`` —
        suitable for unpacking directly into :func:`radarplot` (``**prep``).
    """
    if group_field is None or group_field == '':
        _groups = [None]
    else:
        if groups is None:
            _groups = sorted(data[group_field].dropna().unique())
        else:
            _groups = list(groups)

    n_groups = len(_groups)

    if quantiles is not None:
        n_q = len(quantiles)
        vals = np.full((n_groups, len(fields), n_q), np.nan)
        for i, g in enumerate(_groups):
            subset = data if g is None else data[data[group_field] == g]
            for j, q in enumerate(quantiles):
                vals[i, :, j] = subset[fields].quantile(q).to_numpy()
    else:
        vals = np.full((n_groups, len(fields)), np.nan)
        for i, g in enumerate(_groups):
            subset = data if g is None else data[data[group_field] == g]
            vals[i, :] = subset[fields].mean().to_numpy()

    return dict(vals=vals, fields=fields, groups=_groups, quantiles=quantiles)


# ---------------------------------------------------------------------------
# Parameter resolution helpers
# ---------------------------------------------------------------------------

def _resolve_tick_fields(tick_fields, n_fields):
    """Return a list of field indices that should receive tick labels."""
    if tick_fields == 'all':
        return list(range(n_fields))
    if tick_fields is None:
        return []
    if isinstance(tick_fields, (int, np.integer)):
        return [int(tick_fields)]
    return [int(f) for f in tick_fields]


def _resolve_log(log, fields):
    """Return a boolean mask of length n_fields for log-scaled axes."""
    n = len(fields)
    if isinstance(log, (bool, np.bool_)):
        return np.full(n, bool(log), dtype=bool)
    mask = np.zeros(n, dtype=bool)
    for item in log:
        if isinstance(item, str):
            mask[list(fields).index(item)] = True
        elif isinstance(item, (bool, np.bool_)):
            mask[:] = bool(item)
        else:
            mask[int(item)] = True
    return mask


def _resolve_axis_lim(axis_lim, fields, auto_min, auto_max):
    """Merge user-supplied limits with auto-computed ones.

    *axis_lim* can be:
    - ``None`` — use auto for all fields
    - ``dict`` mapping field name or index to ``(lo, hi)`` — partial override
    - array-like of shape ``(2, n_fields)`` — full explicit limits
    """
    out_min = auto_min.copy()
    out_max = auto_max.copy()

    if axis_lim is None:
        pass
    elif isinstance(axis_lim, dict):
        for k, (lo, hi) in axis_lim.items():
            j = list(fields).index(k) if isinstance(k, str) else int(k)
            if lo is not None:
                out_min[j] = lo
            if hi is not None:
                out_max[j] = hi
    else:
        arr = np.asarray(axis_lim, dtype=float)
        out_min = arr[0].copy()
        out_max = arr[1].copy()

    # Fix zero-span axes
    flat = out_min == out_max
    out_min[flat] -= 0.1 * np.abs(out_min[flat]) + 1e-9
    out_max[flat] += 0.1 * np.abs(out_max[flat]) + 1e-9

    return out_min, out_max


def _resolve_axis_labels(axis_labels, fields):
    """Return display label list, applying any overrides from *axis_labels*.

    *axis_labels* can be:
    - ``None`` — use *fields* as-is
    - ``list`` of str — full replacement (must match length of *fields*)
    - ``dict`` mapping field name or index to display string — partial override
    """
    labels = list(fields)
    if axis_labels is None:
        return labels
    if isinstance(axis_labels, dict):
        for k, v in axis_labels.items():
            j = labels.index(k) if isinstance(k, str) else int(k)
            labels[j] = v
        return labels
    return list(axis_labels)


# ---------------------------------------------------------------------------
# Internal geometry helpers
# ---------------------------------------------------------------------------

def _spoke_angles(n: int, start: float = np.pi / 2, clockwise: bool = True) -> np.ndarray:
    sign = -1.0 if clockwise else 1.0
    return start + sign * np.arange(n) * 2 * np.pi / n


def _to_xy(r, theta):
    return r * np.cos(theta), r * np.sin(theta)


def _normalize_vals(vals, axis_min, axis_max, r_inner, r_outer, log_mask):
    """Map each field's values to [r_inner, r_outer], with per-field log support.

    *vals* shape: (..., n_fields).
    """
    vals = np.asarray(vals, dtype=float)
    result = np.empty_like(vals)
    for j in range(len(axis_min)):
        v = vals[..., j]
        lo, hi = float(axis_min[j]), float(axis_max[j])
        if log_mask[j]:
            lo_s = lo if lo > 0 else 1e-300
            hi_s = max(hi, lo_s * (1 + 1e-9))
            v_n = (np.log10(np.maximum(v, lo_s)) - np.log10(lo_s)) / (
                np.log10(hi_s) - np.log10(lo_s)
            )
        else:
            span = hi - lo if hi != lo else 1.0
            v_n = (v - lo) / span
        result[..., j] = r_inner + np.clip(v_n, 0.0, 1.0) * (r_outer - r_inner)
    return result


def _field_tick_vals(j, axis_min, axis_max, n_ticks, log_mask):
    """Return n_ticks+1 data values for ring tick labels on field *j*."""
    lo, hi = float(axis_min[j]), float(axis_max[j])
    if log_mask[j]:
        lo_s = lo if lo > 0 else 1e-300
        hi_s = max(hi, lo_s * (1 + 1e-9))
        return np.logspace(np.log10(lo_s), np.log10(hi_s), n_ticks + 1)
    return np.linspace(lo, hi, n_ticks + 1)


def _label_alignment(theta: float) -> dict:
    """ha/va alignment pointing away from the centre along spoke *theta*."""
    x, y = np.cos(theta), np.sin(theta)
    ha = 'center' if abs(x) < 0.15 else ('left' if x > 0 else 'right')
    va = 'center' if abs(y) < 0.15 else ('bottom' if y > 0 else 'top')
    return {'ha': ha, 'va': va}


def _draw_web(ax, thetas, n_ticks, r_outer, r_inner, grey):
    """Draw dashed radial spokes and dotted concentric isocurve rings."""
    for t in thetas:
        x0, y0 = _to_xy(r_inner, t)
        x1, y1 = _to_xy(r_outer, t)
        ax.plot([x0, x1], [y0, y1], '--', color=grey, lw=0.8, zorder=0)

    ring_r = np.linspace(r_inner, r_outer, n_ticks + 1)
    theta_ring = np.append(thetas, thetas[0])
    for r in ring_r:
        xr, yr = _to_xy(r, theta_ring)
        ax.plot(xr, yr, ':', color=grey, lw=0.7, zorder=0)

    return ring_r


def _draw_tick_labels(ax, thetas, ring_r, axis_min, axis_max, n_ticks,
                      precision, tick_field_indices, log_mask):
    """Place numeric scale labels along each spoke listed in *tick_field_indices*."""
    pad = 0.025
    for j in tick_field_indices:
        t = thetas[j]
        cx, cy = np.cos(t), np.sin(t)
        tick_vals = _field_tick_vals(j, axis_min, axis_max, n_ticks, log_mask)
        kw = _label_alignment(t)
        for r, tv in zip(ring_r, tick_vals):
            ax.text((r + pad) * cx, (r + pad) * cy,
                    f'{tv:.{precision}g}',
                    fontsize=6.5, color='0.35', zorder=3, **kw)


def _draw_field_labels(ax, thetas, r_outer, labels, pad=0.12):
    """Place bold field-name labels at spoke tips."""
    for t, lbl in zip(thetas, labels):
        x, y = _to_xy(r_outer + pad, t)
        kw = _label_alignment(t)
        ax.text(x, y, lbl, fontsize=9, fontweight='bold', zorder=4, **kw)


# ---------------------------------------------------------------------------
# Public plot function
# ---------------------------------------------------------------------------

_DEFAULT_COLORS = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
    '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
    '#bcbd22', '#17becf',
]


def radarplot(
    vals,
    fields: list[str],
    groups=None,
    quantiles: list[float] | None = None,
    style: str = 'web',
    ax: plt.Axes | None = None,
    # axis configuration
    n_ticks: int = 4,
    precision: int = 3,
    axis_lim=None,
    log=False,
    axis_labels=None,
    tick_fields='all',
    # web-style options
    colors=None,
    alpha_fill: float = 0.15,
    clockwise: bool = True,
    legend_loc: str = 'upper right',
    linewidth: float = 1.0,
    # series-style options
    cmap: str = 'viridis',
) -> plt.Axes:
    """Draw a radar (spider) or sequential radar plot.

    Parameters
    ----------
    vals : array-like
        For ``style='web'``: shape ``(n_groups, n_fields)`` or
        ``(n_groups, n_fields, n_quantiles)``.
        For ``style='series'``: shape ``(n_obs, n_fields)`` — each row is
        one lap around the spokes.
    fields : list of str
        Data field names — also used as spoke labels unless *axis_labels*
        overrides them.
    groups : list, optional
        Group names for the legend (``style='web'``).
    quantiles : list of float, optional
        Quantile levels matching the third axis of *vals*.
    style : ``'web'`` | ``'series'``
        Plot style.  Default ``'web'``.
    ax : matplotlib.axes.Axes, optional
        Target axes.  Created if *None*.
    n_ticks : int
        Number of concentric ring intervals.  Default 4.
    precision : int
        Significant figures for tick labels.  Default 3.
    axis_lim : None | dict | array-like (2, n_fields)
        Per-axis data limits.

        - ``None`` — auto-compute from data (default).
        - ``dict`` mapping field name or index to ``(lo, hi)`` — partial
          override; unspecified fields remain auto-computed.  Pass ``None``
          for either bound to keep it auto, e.g.
          ``{'SiO2': (40, 80), 'TiO2': (None, 5)}``.
        - Array-like of shape ``(2, n_fields)`` — full explicit ``[min; max]``.
    log : bool | list of str | list of int
        Log-scale axes.

        - ``False`` (default) — all linear.
        - ``True`` — all log.
        - List of field names or indices — only those fields use log scale.
    axis_labels : None | list of str | dict
        Display labels for spokes, overriding *fields*.

        - ``None`` — use *fields* as-is (default).
        - ``list`` — full replacement, same length as *fields*.
        - ``dict`` mapping field name or index to display string — partial
          override; e.g. ``{'TiO2': 'TiO₂ (log)'}``.
    tick_fields : ``'all'`` | None | int | list of int
        Which spokes get numeric tick labels.

        - ``'all'`` — label every spoke (default).
        - ``None`` — no tick labels.
        - ``int`` or list of int — label only those spoke indices.
    colors : list, optional
        Per-group colours for ``style='web'``.  Cycles through the default
        Matplotlib colour sequence if *None*.
    alpha_fill : float
        Transparency of quantile-band fills.  Default 0.15.
    clockwise : bool
        Spoke direction from the top.  Default ``True``.
    legend_loc : str
        Legend location.  Default ``'upper right'``.
    linewidth : float
        Data line width.  Default 1.0.
    cmap : str
        Colormap for ``style='series'``.  Default ``'viridis'``.

    Returns
    -------
    matplotlib.axes.Axes
    """
    vals = np.asarray(vals, dtype=float)
    style = style.lower()
    n_fields = len(fields)

    log_mask = _resolve_log(log, fields)
    tick_idx = _resolve_tick_fields(tick_fields, n_fields)
    display_labels = _resolve_axis_labels(axis_labels, fields)

    if ax is None:
        _, ax = plt.subplots(figsize=(7, 7))

    if style == 'web':
        return _radarplot_web(
            ax, vals, fields, display_labels, groups, quantiles,
            n_ticks, precision, axis_lim, log_mask, tick_idx,
            colors, alpha_fill, clockwise, legend_loc, linewidth,
        )
    elif style == 'series':
        return _radarplot_series(
            ax, vals, fields, display_labels,
            n_ticks, precision, axis_lim, log_mask, tick_idx,
            clockwise, cmap, linewidth,
        )
    else:
        raise ValueError(f"Unknown style '{style}'. Choose 'web' or 'series'.")


# ---------------------------------------------------------------------------
# Web style
# ---------------------------------------------------------------------------

def _radarplot_web(ax, vals, fields, display_labels, groups, quantiles,
                   n_ticks, precision, axis_lim, log_mask, tick_idx,
                   colors, alpha_fill, clockwise, legend_loc, linewidth):
    n_fields = len(fields)
    grey = '0.6'

    # Normalise to (n_groups, n_fields, n_q)
    if vals.ndim == 2:
        vals3 = vals[:, :, np.newaxis]
        n_groups = vals.shape[0]
        n_q = 1
    else:
        vals3 = vals
        n_groups, _, n_q = vals.shape

    # Auto-compute limits from the full quantile range, then apply overrides
    auto_min = np.nanmin(vals3[:, :, 0], axis=0)
    auto_max = np.nanmax(vals3[:, :, -1], axis=0)
    axis_min, axis_max = _resolve_axis_lim(axis_lim, fields, auto_min, auto_max)

    r_inner = 1.0 / (n_ticks + 1)
    r_outer = 1.0

    thetas = _spoke_angles(n_fields, clockwise=clockwise)
    ring_r = _draw_web(ax, thetas, n_ticks, r_outer, r_inner, grey)
    _draw_field_labels(ax, thetas, r_outer, display_labels)
    _draw_tick_labels(ax, thetas, ring_r, axis_min, axis_max,
                      n_ticks, precision, tick_idx, log_mask)

    if colors is None:
        colors = (_DEFAULT_COLORS * ((n_groups // len(_DEFAULT_COLORS)) + 1))[:n_groups]

    theta_closed = np.append(thetas, thetas[0])
    legend_handles = []

    for i in range(n_groups):
        c = colors[i]
        _draw_group(ax, vals3[i], axis_min, axis_max, r_inner, r_outer,
                    log_mask, n_q, theta_closed, c, alpha_fill, linewidth)
        label = str(groups[i]) if (groups is not None and groups[i] is not None) else None
        legend_handles.append(mpatches.Patch(color=c, label=label, alpha=0.7))

    if groups is not None and any(g is not None for g in groups):
        ax.legend(handles=legend_handles, loc=legend_loc, frameon=False, fontsize=8)

    _finish_ax(ax, r_outer)
    return ax


def _draw_group(ax, g_vals, axis_min, axis_max, r_inner, r_outer,
                log_mask, n_q, theta_closed, color, alpha_fill, linewidth):
    """Draw one group's polygon(s) for the web style."""

    def _r(q_idx):
        return _normalize_vals(
            g_vals[:, q_idx], axis_min, axis_max, r_inner, r_outer, log_mask
        )

    def _line(r, **kw):
        rc = np.append(r, r[0])
        x, y = rc * np.cos(theta_closed), rc * np.sin(theta_closed)
        ax.plot(x, y, **kw)

    if n_q == 1:
        r0 = _r(0)
        _line(r0, color=color, lw=linewidth, zorder=2)
        ax.fill(*(np.append(r0, r0[0]) * np.array([np.cos(theta_closed),
                                                     np.sin(theta_closed)])),
                color=color, alpha=alpha_fill, zorder=1)

    elif n_q == 2:
        _fill_band(ax, theta_closed, _r(0), _r(1), color, alpha_fill)
        _line(_r(1), color=color, lw=linewidth, zorder=2)

    elif n_q == 3:
        _fill_band(ax, theta_closed, _r(0), _r(2), color, alpha_fill)
        _line(_r(1), color=color, lw=linewidth, zorder=2)

    elif n_q == 5:
        for ri in (_r(0), _r(4)):
            _line(ri, color=color, lw=linewidth * 0.5, ls='--', alpha=0.5, zorder=2)
        _fill_band(ax, theta_closed, _r(1), _r(3), color, alpha_fill)
        _line(_r(2), color=color, lw=linewidth, zorder=2)

    else:
        r_all = [_r(k) for k in range(n_q)]
        _fill_band(ax, theta_closed, r_all[0], r_all[-1], color, alpha_fill)
        for ri in r_all:
            _line(ri, color=color, lw=linewidth * 0.5, zorder=2)


def _fill_band(ax, theta_closed, r_lo, r_hi, color, alpha):
    """Fill the annular region between two radius arrays."""
    def _xy(r):
        rc = np.append(r, r[0])
        return rc * np.cos(theta_closed), rc * np.sin(theta_closed)

    xo, yo = _xy(r_hi)
    xi, yi = _xy(r_lo)
    x = np.concatenate([xo, xi[::-1], [xo[0]]])
    y = np.concatenate([yo, yi[::-1], [yo[0]]])
    ax.fill(x, y, color=color, alpha=alpha, zorder=1, ec='none')


# ---------------------------------------------------------------------------
# Series style
# ---------------------------------------------------------------------------

def _radarplot_series(ax, vals, fields, display_labels,
                      n_ticks, precision, axis_lim, log_mask, tick_idx,
                      clockwise, cmap, linewidth):
    """Sequential radar: each row of *vals* is one lap; consecutive laps connect."""
    n_obs, n_fields = vals.shape
    grey = '0.6'

    auto_min = np.nanmin(vals, axis=0)
    auto_max = np.nanmax(vals, axis=0)
    axis_min, axis_max = _resolve_axis_lim(axis_lim, fields, auto_min, auto_max)

    r_inner = 1.0 / (n_ticks + 1)
    r_outer = 1.0

    thetas = _spoke_angles(n_fields, clockwise=clockwise)
    ring_r = _draw_web(ax, thetas, n_ticks, r_outer, r_inner, grey)
    _draw_field_labels(ax, thetas, r_outer, display_labels)
    _draw_tick_labels(ax, thetas, ring_r, axis_min, axis_max,
                      n_ticks, precision, tick_idx, log_mask)

    r_norm = _normalize_vals(vals, axis_min, axis_max, r_inner, r_outer, log_mask)

    pts_x, pts_y = [], []
    for k in range(n_obs):
        for j in range(n_fields):
            pts_x.append(r_norm[k, j] * np.cos(thetas[j]))
            pts_y.append(r_norm[k, j] * np.sin(thetas[j]))
        # Close each lap back to the first spoke
        pts_x.append(r_norm[k, 0] * np.cos(thetas[0]))
        pts_y.append(r_norm[k, 0] * np.sin(thetas[0]))

    pts = np.column_stack([pts_x, pts_y])
    segments = np.stack([pts[:-1], pts[1:]], axis=1)
    color_vals = np.linspace(0.0, 1.0, len(segments))

    lc = LineCollection(segments, cmap=cmap, norm=mcolors.Normalize(0, 1),
                        linewidth=linewidth, zorder=2, alpha=0.85)
    lc.set_array(color_vals)
    ax.add_collection(lc)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=mcolors.Normalize(0, n_obs - 1))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.6, pad=0.04, aspect=20)
    cbar.set_label('Observation', fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    _finish_ax(ax, r_outer)
    return ax


# ---------------------------------------------------------------------------
# Shared axis finish
# ---------------------------------------------------------------------------

def _finish_ax(ax, r_outer):
    lim = r_outer * 1.35
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect('equal')
    ax.axis('off')
