"""
Heat-production vs seismic-velocity plotting utilities.

Ports of matlab/plotting/hpVplot.m, hpax.m, and hpVhist.m.
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from scipy.odr import ODR, Model, RealData

from src.plotting.whisker import whisker


# ── axis / field configuration ────────────────────────────────────────────────

_FIELD_DEFAULTS = {
    'p_velocity': dict(
        xlabel='Estimated P-Velocity [km s⁻¹]',
        xlim=(5.8, 8.2),
        xfit=(6.0, 7.4),
        xshift=6.0,
        bin_edges=np.arange(5.8, 8.21, 0.1),
    ),
    's_velocity': dict(
        xlabel='Estimated S-Velocity [km s⁻¹]',
        xlim=(3.65, 4.6),
        xfit=(3.65, 4.05),
        xshift=3.65,
        bin_edges=np.arange(3.6, 4.61, 0.05),
    ),
    'vmoho': dict(
        xlabel='Crustal Thickness [km]',
        xlim=(0.0, 75.0),
        xfit=(0.0, 75.0),
        xshift=0.0,
        bin_edges=np.arange(0.0, 76.0, 5.0),
    ),
}


# ── orthogonal distance regression (errors in both x and y) ──────────────────

def _odr_fit(x, y, sigma_x, sigma_y):
    """Linear ODR fit: y = m[0]*x + m[1] with errors on both axes.

    Matches the 4-argument ``wllsq`` called inside ``hpVplot.m``.

    Parameters
    ----------
    x, y : array-like, 1-D
    sigma_x, sigma_y : array-like
        1-sigma uncertainties on x and y respectively.

    Returns
    -------
    m : ndarray, shape (2,)
        [slope, intercept]
    m_err : ndarray, shape (2,)
        1-sigma uncertainties on m.
    rsq : float
        R² from the ODR residuals.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    sx = np.asarray(sigma_x, dtype=float)
    sy = np.asarray(sigma_y, dtype=float)

    # guard against zero / NaN sigmas
    sx = np.where((sx > 0) & np.isfinite(sx), sx, 0.15 * np.abs(x))
    sy = np.where((sy > 0) & np.isfinite(sy), sy, 0.15 * np.abs(y))
    sx[sx == 0] = 1e-6
    sy[sy == 0] = 1e-6

    mask = np.isfinite(x) & np.isfinite(y)
    x, y, sx, sy = x[mask], y[mask], sx[mask], sy[mask]

    def linear(beta, t):
        return beta[0] * t + beta[1]

    # seed from ordinary least squares
    beta0 = np.polyfit(x, y, 1)
    odr_data  = RealData(x, y, sx=sx, sy=sy)
    odr_model = Model(linear)
    out = ODR(odr_data, odr_model, beta0=beta0).run()

    m     = out.beta           # [slope, intercept]
    m_err = out.sd_beta

    ss_res = np.sum(((y - linear(m, x)) / sy)**2)
    ss_tot = np.sum(((y - np.mean(y)) / sy)**2)
    rsq = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

    return m, m_err, rsq


# ── axis formatter ────────────────────────────────────────────────────────────

def hp_ax(lim=(-2, 1), axis='y', ax=None):
    """Configure an axis as a log-scale heat-production axis.

    Mirrors MATLAB ``hpax``.  The axis displays values in μW m⁻³ but the
    underlying data is log₁₀(heat production).

    Parameters
    ----------
    lim : (int, int)
        (min, max) exponents, e.g. ``(-2, 1)`` gives 0.01 – 10 μW m⁻³.
    axis : {'x', 'y'}
    ax : matplotlib Axes, optional
    """
    if ax is None:
        ax = plt.gca()

    minor_exp = np.log10(np.arange(2, 10))
    ticks, labels = [], []
    for exp in range(lim[0], lim[1] + 1):
        ticks.append(exp)
        labels.append(f'10^{exp}' if exp != 0 else '1')
        for me in minor_exp:
            ticks.append(exp + me)
            labels.append('')

    hp_label = 'Heat Production [μW m⁻³]'
    if axis == 'x':
        ax.set_xlabel(hp_label)
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels)
        ax.set_xlim(lim)
    else:
        ax.set_ylabel(hp_label)
        ax.set_yticks(ticks)
        ax.set_yticklabels(labels)
        ax.set_ylim(lim)
    ax.set_box_aspect(None)


# ── main plot function ────────────────────────────────────────────────────────

def hp_velocity_plot(
    V,
    A,
    field='p_velocity',
    color=None,
    bin_edges=None,
    min_count=9,
    ylim=(-2, 1),
    ax=None,
    A_raw=None,
):
    """Plot log₁₀(heat production) vs seismic velocity with binned quantiles.

    Mirrors MATLAB ``hpVplot``.  Each velocity bin is rendered as a
    box-and-whisker (2.5 / 25 / 50 / 75 / 97.5 percentiles).  A weighted
    least-squares line is fitted to the bin medians within *xfit*.

    Parameters
    ----------
    V : array-like
        Velocity values (km s⁻¹ for P/S, km for Moho depth).
    A : array-like
        Heat production values (μW m⁻³, positive).
    field : {'p_velocity', 's_velocity', 'vmoho'}
    color : colour spec, optional
        Default is a steel blue ``(54, 169, 225) / 255``.
    bin_edges : array-like, optional
        Override the default bin edges for *field*.
    min_count : int
        Minimum number of points in a bin to include in the fit.
    ylim : (float, float)
        log₁₀ y-axis limits.
    ax : matplotlib Axes, optional
    A_raw : array-like, optional
        Raw heat-production values with the BDL convention (positive =
        detected, negative = BDL with known detection limit, zero = BDL
        with unknown detection limit).  When provided, per-bin Y quantiles
        are estimated via ``gausscensor`` rather than ``numpy.quantile``.
        Must have the same length as *A*.

    Returns
    -------
    m : ndarray, shape (2,)
        [slope, intercept] of weighted fit.
    m_err : ndarray, shape (2,)
        1-sigma uncertainties on m.
    C : float
        Signed correlation coefficient (sign taken from slope).
    bin_edges : ndarray
        The bin edges actually used (for histogram reuse).
    """
    if color is None:
        color = np.array([54, 169, 225]) / 255.0

    cfg = _FIELD_DEFAULTS[field]
    if bin_edges is None:
        bin_edges = cfg['bin_edges'].copy()
    bin_edges = np.asarray(bin_edges, dtype=float)

    xfit   = cfg['xfit']
    xshift = cfg['xshift']

    V = np.asarray(V, dtype=float).ravel()
    A = np.asarray(A, dtype=float).ravel()
    A_raw_arr = np.asarray(A_raw, dtype=float).ravel() if A_raw is not None else None

    if ax is None:
        ax = plt.gca()

    # ── bin data into cell-array style lists ──────────────────────────────────
    n_bins = len(bin_edges) - 1
    V_bins, A_bins, A_raw_bins = [], [], []
    N = np.zeros(n_bins, dtype=int)
    for i in range(n_bins):
        mask = (bin_edges[i] <= V) & (V < bin_edges[i + 1]) & (A > 0)
        N[i] = mask.sum()
        V_bins.append(V[mask])
        A_bins.append(A[mask])
        if A_raw_arr is not None:
            A_raw_bins.append(A_raw_arr[mask])

    # ── draw 2-D box-and-whisker plot (whisker.m behaviour) ──────────────────
    X, Y = whisker(
        V_bins, A_bins,
        scale='log',
        color=color,
        plot=True,
        ax=ax,
        y_raw=A_raw_bins if A_raw_arr is not None else None,
    )

    # ── ODR fit to bin medians (errors in both x and y) ──────────────────────
    # 1-SE estimated from IQR / (1.35 * sqrt(N))  (robust, matches hpVplot.m)
    iqr_x = X[:, 3] - X[:, 1]
    iqr_y = Y[:, 3] - Y[:, 1]
    with np.errstate(invalid='ignore', divide='ignore'):
        sqrtN    = np.sqrt(N.astype(float))
        sigma_x  = iqr_x / (1.35 * sqrtN)
        sigma_y  = iqr_y / (1.35 * sqrtN)

    fit_mask = (
        (xfit[0] <= X[:, 2]) & (X[:, 2] <= xfit[1])
        & np.isfinite(Y[:, 2])
        & np.isfinite(sigma_x) & np.isfinite(sigma_y)
        & (N >= min_count)
    )

    nan2 = np.array([np.nan, np.nan])
    if fit_mask.sum() >= 3:
        m, m_err, rsq = _odr_fit(
            X[fit_mask, 2] - xshift,
            Y[fit_mask, 2],
            sigma_x[fit_mask],
            sigma_y[fit_mask],
        )
        C = np.sign(m[0]) * np.sqrt(max(rsq, 0.0))

        Xplot = np.linspace(bin_edges[0], bin_edges[-1], 200)
        ax.plot(Xplot, m[0] * (Xplot - xshift) + m[1], '-', color=color)
    else:
        m, m_err, C = nan2, nan2, np.nan

    hp_ax(ylim, axis='y', ax=ax)
    ax.set_xlabel(cfg['xlabel'])
    ax.set_xlim(cfg['xlim'])
    ax.set_box_aspect(None)

    return m, m_err, C, bin_edges


# ── histogram panel ───────────────────────────────────────────────────────────

def hp_velocity_hist(sio2, V, A, rho, sio2_edges, axes=None):
    """SiO₂-binned histograms of log₁₀(A), density, and velocity.

    Mirrors MATLAB ``hpVhist``.  Produces a grid of ``n_bins × 3``
    subplots: log₁₀(heat production) | density | P-velocity, one row per
    SiO₂ interval (bottom row = lowest SiO₂).

    Parameters
    ----------
    sio2 : array-like
        SiO₂ wt% for each sample.
    V : array-like
        P-velocity (km s⁻¹).
    A : array-like
        Heat production (μW m⁻³, positive).
    rho : array-like
        Density (kg m⁻³).
    sio2_edges : array-like
        SiO₂ bin edges, e.g. ``[45, 52, 63, 75]``.
    axes : array of Axes, shape (n_bins, 3), optional
        Supply pre-created axes; otherwise a new figure is created.

    Returns
    -------
    fig : matplotlib Figure
    axes : ndarray of Axes
    """
    sio2 = np.asarray(sio2, dtype=float).ravel()
    V    = np.asarray(V,    dtype=float).ravel()
    A    = np.asarray(A,    dtype=float).ravel()
    rho  = np.asarray(rho,  dtype=float).ravel()
    sio2_edges = np.asarray(sio2_edges, dtype=float)

    nr = len(sio2_edges) - 1

    if axes is None:
        fig, axes = plt.subplots(nr, 3, figsize=(9, 3 * nr))
        axes = np.atleast_2d(axes)
    else:
        fig = axes.flat[0].get_figure()

    logA = np.where(A > 0, np.log10(A), np.nan)

    for i in range(nr):
        row = nr - 1 - i          # bottom row = lowest SiO₂ (matches MATLAB)
        mask = (sio2_edges[i] <= sio2) & (sio2 < sio2_edges[i + 1])

        # log₁₀(A)
        ax = axes[row, 0]
        ax.hist(logA[mask], bins=np.arange(-3, 3.1, 0.25),
                density=True, histtype='step')
        ax.set_xlim(-2, 2)
        if i == 0:
            ax.set_xlabel('Heat Production [μW m⁻³]')
            ax.set_xticks(range(-3, 4))
            ax.set_xticklabels([f'10^{e}' if e != 0 else '1' for e in range(-3, 4)])
        elif i == nr - 1:
            ax.xaxis.set_label_position('top')
            ax.xaxis.tick_top()
        else:
            ax.axis('off')

        # density
        ax = axes[row, 1]
        ax.hist(rho[mask], bins=np.arange(2500, 3501, 50),
                density=True, histtype='step')
        ax.set_xlim(2500, 3500)
        if i == 0:
            ax.set_xlabel('Density [kg m⁻³]')
            ax.set_xticks(range(2600, 3401, 200))
        elif i == nr - 1:
            ax.xaxis.set_label_position('top')
            ax.xaxis.tick_top()
            ax.set_xticks(range(2600, 3401, 200))
        else:
            ax.axis('off')

        # P-velocity
        ax = axes[row, 2]
        ax.hist(V[mask], bins=np.arange(5.8, 8.41, 0.1),
                density=True, histtype='step')
        ax.set_xlim(5.8, 8.4)
        if i == 0:
            ax.set_xlabel('Estimated P-Velocity [km s⁻¹]')
        elif i == nr - 1:
            ax.xaxis.set_label_position('top')
            ax.xaxis.tick_top()
        else:
            ax.axis('off')

    fig.tight_layout()
    return fig, axes
