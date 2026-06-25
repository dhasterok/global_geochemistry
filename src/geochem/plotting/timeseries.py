"""
Binned time series ribbon plot for geochemical data through geologic time.

``age_series(data, el1, el2=None, ...)``
    Bins data by age and draws a step-function ribbon showing the median and
    interquartile range (IQR), with optional 2.5–97.5 % extent lines.
    Quantiles are estimated via :func:`~.processing.censored.gausscensor` to
    account for below-detection-limit values.

    Covers the functionality of both MATLAB ``elratio_vs_age.m`` and the
    plotting section of ``decay_correction.m``.

Typical use
-----------
Single element vs. age::

    ax = age_series(data, 'th_ppm', scale='log')

Element ratio vs. age with decay correction::

    ax = age_series(data, 'th_ppm', el2='u_ppm',
                    decay_correct=True, scale='log')

Overlay two ribbons on the same axes::

    fig, ax = plt.subplots()
    age_series(data_ign, 'u_ppm',  ax=ax, color='steelblue', label='Igneous')
    age_series(data_sed, 'u_ppm',  ax=ax, color='tomato',    label='Sedimentary')
    ax.legend()
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

from global_geochemistry.geochem.processing.censored import gausscensor
from global_geochemistry.geochem.processing.decay    import decaycorrect, DECAY_ELEMENTS


# Default age bin edges (Ma): 0 → 4000 in 500 Ma steps
_DEFAULT_EDGES = np.arange(0, 4001, 500, dtype=float)

_QUANTS = [0.025, 0.25, 0.5, 0.75, 0.975]


def age_series(
    data: pd.DataFrame,
    el1: str,
    el2: str | None = None,
    bin_edges=None,
    age_col: str = 'avg_age',
    scale: str = 'log',
    decay_correct: bool = False,
    decay_ref: str = 'r88',
    color='steelblue',
    alpha_band: float = 0.4,
    min_n: int = 10,
    label: str | None = None,
    ax=None,
) -> tuple[plt.Axes, dict]:
    """Plot element (or ratio) concentrations binned by age.

    Parameters
    ----------
    data : pandas.DataFrame
        Must contain *el1* (and *el2* if provided) and *age_col*.
    el1 : str
        Column name of the numerator element (or single element if *el2* is
        None).
    el2 : str, optional
        Column name of the denominator element.  When provided, plots
        ``el1 / el2``.
    bin_edges : array-like, optional
        Age bin edges in Ma.  Defaults to ``0, 500, 1000, …, 4000 Ma``.
    age_col : str
        Name of the age column (default ``'avg_age'``).
    scale : {'log', 'linear'}
        Y-axis scale.  ``'log'`` uses :func:`~.processing.censored.gausscensor`
        in log₁₀ space so that quantile uncertainty is symmetric on a log axis.
    decay_correct : bool
        If ``True``, apply radioactive decay correction to each sample before
        binning.  Only meaningful for K, Rb, Sm, Th, U.
    decay_ref : str
        Reference for decay constants forwarded to
        :func:`~.processing.decay.decaycorrect` (default ``'r88'``).
    color : matplotlib color
        Ribbon fill and line colour.
    alpha_band : float
        Opacity of the IQR fill band (default 0.4).
    min_n : int
        Minimum samples per bin; bins with fewer samples are skipped (default 10).
    label : str, optional
        Legend label for this series.
    ax : matplotlib.axes.Axes, optional
        Target axes.  A new figure is created if *None*.

    Returns
    -------
    ax : matplotlib.axes.Axes
    stats : dict
        ``{'bin_edges': ndarray, 'xq': (5, N_bins), 'yq': (5, N_bins), 'n': (N_bins,)}``.
        ``xq`` is the age quantiles; ``yq`` is the element/ratio quantiles.
        Both are in original (linear) units regardless of *scale*.

    Examples
    --------
    >>> ax, stats = age_series(data, 'th_ppm', scale='log')
    >>> ax, stats = age_series(data, 'th_ppm', 'u_ppm',
    ...                        decay_correct=True, color='tomato')
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 4))

    edges = np.asarray(bin_edges if bin_edges is not None else _DEFAULT_EDGES,
                       dtype=float)
    n_bins = len(edges) - 1

    # --- extract and filter data ---
    sub = data.copy()
    missing_cols = [c for c in [el1, age_col] + ([el2] if el2 else [])
                    if c not in sub.columns]
    if missing_cols:
        raise ValueError(f"Columns not found in data: {missing_cols}")

    v1 = sub[el1].to_numpy(float)
    v2 = sub[el2].to_numpy(float) if el2 else None
    age = sub[age_col].to_numpy(float)

    if v2 is not None:
        valid = (v1 > 0) & (v2 > 0) & np.isfinite(age)
    else:
        valid = (v1 > 0) & np.isfinite(age)
    v1, age_v = v1[valid], age[valid]
    if v2 is not None:
        v2 = v2[valid]

    # --- optional decay correction ---
    if decay_correct:
        el1_key = el1.lower().removesuffix('_ppm')
        if el1_key in DECAY_ELEMENTS:
            v1 = decaycorrect(el1, v1, age_v, ref=decay_ref)
        if v2 is not None:
            el2_key = el2.lower().removesuffix('_ppm')
            if el2_key in DECAY_ELEMENTS:
                v2 = decaycorrect(el2, v2, age_v, ref=decay_ref)

    # Compute values to bin
    values = v1 / v2 if v2 is not None else v1

    # --- binning and quantile estimation ---
    censor_scale = 'log10' if scale == 'log' else 'linear'
    xq_list, yq_list, n_list = [], [], []

    for i in range(n_bins):
        mask = (edges[i] <= age_v) & (age_v < edges[i + 1])
        n = mask.sum()
        n_list.append(n)
        if n < min_n:
            xq_list.append(np.full(5, np.nan))
            yq_list.append(np.full(5, np.nan))
            continue

        age_bin = age_v[mask]
        val_bin = values[mask]

        # Age quantiles (always linear, no BDL handling needed)
        xq_list.append(np.quantile(age_bin, _QUANTS))

        # Value quantiles via gausscensor
        result = gausscensor(val_bin, scale=censor_scale, q=_QUANTS)
        if not isinstance(result, tuple):
            yq_list.append(np.full(5, np.nan))
            continue
        _, qvals = result
        q = qvals[:, 1]   # empirical quantiles in transformed space
        if scale == 'log':
            q = 10.0 ** q  # back-transform to original units
        yq_list.append(q)

    xq = np.column_stack(xq_list)   # (5, N_bins)
    yq = np.column_stack(yq_list)
    n  = np.array(n_list)

    # --- draw ribbon ---
    _draw_ribbon(ax, edges, yq, color=color, alpha_band=alpha_band, label=label)

    # --- axis formatting ---
    ax.set_xscale('linear')
    ax.set_yscale(scale)
    ax.set_xlim(edges[0], edges[-1])
    _set_age_ticks(ax, edges)
    ax.set_xlabel('Age (Ga)')
    ax.set_ylabel(_make_ylabel(el1, el2))
    ax.set_box_aspect(None)

    stats = {'bin_edges': edges, 'xq': xq, 'yq': yq, 'n': n}
    return ax, stats


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _draw_ribbon(ax, edges, yq, color, alpha_band, label):
    """Draw step-function ribbon: IQR band + median + extent lines.

    yq : (5, N_bins) — rows are [p2.5, p25, p50, p75, p97.5]
    edges : (N_bins+1,) bin edges
    """
    n_bins = yq.shape[1]

    # Build step-function x coordinates: duplicate each edge so we get
    # horizontal segments across each bin
    x_steps = np.empty(2 * n_bins)
    x_steps[0::2] = edges[:-1]
    x_steps[1::2] = edges[1:]

    def _step(q_row):
        """Repeat each bin value twice to match x_steps."""
        return np.repeat(q_row, 2)

    q025  = _step(yq[0])
    q25   = _step(yq[1])
    q50   = _step(yq[2])
    q75   = _step(yq[3])
    q975  = _step(yq[4])

    # IQR fill band (Q25–Q75)
    ax.fill_between(x_steps, q25, q75,
                    color=color, alpha=alpha_band, linewidth=0,
                    label=label)

    # Median line
    ax.plot(x_steps, q50, color=color, linewidth=1.5)

    # 2.5–97.5 % extent lines (thin dashed)
    ax.plot(x_steps, q025, color=color, linewidth=0.75, linestyle='--')
    ax.plot(x_steps, q975, color=color, linewidth=0.75, linestyle='--')


def _set_age_ticks(ax, edges):
    """Label x-axis in Ga with 500 Ma major ticks."""
    x_min, x_max = edges[0], edges[-1]
    tick_ma = np.arange(
        np.ceil(x_min / 500) * 500,
        x_max + 1,
        500,
        dtype=float,
    )
    ax.set_xticks(tick_ma)
    ax.xaxis.set_major_formatter(
        mticker.FuncFormatter(lambda v, _: f'{v / 1000:g}')
    )


def _make_ylabel(el1: str, el2: str | None) -> str:
    """Build a y-axis label from element column name(s)."""
    def _fmt(col):
        c = col.removesuffix('_ppm').removesuffix('_wtper')
        return c.upper() if len(c) <= 2 else c

    if el2 is None:
        return f'{_fmt(el1)} (ppm)'
    return f'{_fmt(el1)}/{_fmt(el2)}'
