"""
Harker variation diagrams.

Plots nine major-oxide panels (TiO2, Al2O3, FeO_T, MgO, CaO, MnO,
Na2O, K2O, P2O5) versus SiO2 in a 3×3 grid.

Two modes are supported:
- ``'scatter'``   — coloured scatter plot (default)
- ``'histogram'`` — 2-D density heatmap (log or linear colour scale)

Reference:
    Harker, A. (1909) The Natural History of Igneous Rocks.
    Methuen & Co., London.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


# Default oxide panels: (column name, axis label, igneous y-max, sed y-max)
_OXIDES = [
    ('tio2',    'TiO$_2$',  5,   3),
    ('al2o3',   'Al$_2$O$_3$', 30, 30),
    ('feo_tot', 'FeO$_T$',  20,  20),
    ('mgo',     'MgO',      20,  20),
    ('cao',     'CaO',      20,  20),
    ('mno',     'MnO',      10,  80),
    ('na2o',    'Na$_2$O',  10,  10),
    ('k2o',     'K$_2$O',   10,  10),
    ('p2o5',    'P$_2$O$_5$', 2.5, 2.5),
]

_FILL_POLY = np.array([[0, 100, 100, 0], [100, 0, 100, 100]])


def harker(
    data,
    plot_type='scatter',
    color=None,
    marker='o',
    alpha=1.0,
    n_bins=25,
    bin_scale='log',
    rock_type='igneous',
    fig=None,
    axes=None,
):
    """Create a 3×3 Harker diagram grid.

    Parameters
    ----------
    data : pandas.DataFrame
        Must contain ``'sio2'`` and at least one of the standard oxide
        columns (see :attr:`_OXIDES`).
    plot_type : {'scatter', 'histogram'}
        Scatter plot (default) or 2-D density heatmap.
    color : array-like, optional
        Scatter only.  A single RGB triplet, an array of length N for
        a sequential colour scale, or an (N, 3) array of per-point colours.
    marker : str
        Scatter marker style (default ``'o'``).
    alpha : float
        Scatter marker opacity 0–1 (default 1).
    n_bins : int
        Number of bins along each axis for histogram mode (default 25).
    bin_scale : {'log', 'linear'}
        Colour scale for histogram mode (default ``'log'``).
    rock_type : {'igneous', 'sedimentary'}
        Controls the default y-axis limits for each panel.
    fig : matplotlib.figure.Figure, optional
        Existing figure to draw into.  If None a new figure is created.
    axes : array-like of Axes, optional
        Nine axes in row-major order.  If provided, *fig* is ignored.

    Returns
    -------
    fig : matplotlib.figure.Figure
    axes : list of matplotlib.axes.Axes (length 9)
    """
    if 'sio2' not in data.columns:
        raise ValueError("data must contain 'sio2'")

    sio2 = data['sio2'].to_numpy(float)
    valid_sio2 = np.isfinite(sio2)
    sio2_range = (
        np.floor(np.nanmin(sio2[valid_sio2]) / 5) * 5,
        np.ceil(np.nanmax(sio2[valid_sio2]) / 5) * 5,
    )

    if axes is None:
        if fig is None:
            fig = plt.figure(figsize=(12, 10))
        gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.45, wspace=0.35)
        axes = [fig.add_subplot(gs[i // 3, i % 3]) for i in range(9)]
    else:
        axes = list(axes)
        fig = axes[0].get_figure()

    use_sed = rock_type == 'sedimentary'

    for ax, (col, ylabel, ymax_ign, ymax_sed) in zip(axes, _OXIDES):
        if col not in data.columns:
            ax.set_visible(False)
            continue

        y = data[col].to_numpy(float)
        ymax_default = ymax_sed if use_sed else ymax_ign
        # Use 130th percentile of the 95th-quantile to set the y limit
        ymax = _auto_ymax(y, ymax_default)
        ylim = (0, ymax)

        ax.set_xlabel('SiO$_2$ (wt%)')
        ax.set_ylabel(f'{ylabel} (wt%)')
        ax.set_xlim(sio2_range)
        ax.set_ylim(ylim)
        ax.tick_params(direction='out')
        ax.set_box_aspect(None)

        if plot_type == 'scatter':
            _scatter_panel(ax, sio2, y, color, marker, alpha)
        elif plot_type == 'histogram':
            _histogram_panel(ax, sio2, y, sio2_range, ylim, n_bins, bin_scale)
        else:
            raise ValueError(f"Unknown plot_type '{plot_type}'")

    return fig, axes


def _auto_ymax(y, default):
    finite = y[np.isfinite(y) & (y >= 0)]
    if len(finite) == 0:
        return default
    q95 = np.quantile(finite, 0.95)
    return max(np.ceil(1.3 * q95), default)


def _scatter_panel(ax, sio2, y, color, marker, alpha):
    """Draw scatter plot with forbidden region fill."""
    # Grey forbidden region (above the total-oxide = 100% line)
    ax.fill(
        [0, 100, 100, 0], [100, 0, 100, 100],
        color=[0.85, 0.85, 0.85], ec='none', zorder=0,
    )

    valid = np.isfinite(sio2) & np.isfinite(y)
    if not valid.any():
        return

    x_v, y_v = sio2[valid], y[valid]
    n = len(x_v)

    if color is None:
        c = [[0.5, 0.5, 0.5]] * n
        show_cbar = False
    else:
        c = np.asarray(color)
        if c.ndim == 1 and c.shape[0] == 3 and n != 3:
            # Single RGB triplet for all points
            c = [c] * n
            show_cbar = False
        elif c.ndim == 1 and len(c) == len(sio2):
            c = c[valid]
            show_cbar = True
        elif c.ndim == 2 and c.shape == (len(sio2), 3):
            c = c[valid]
            show_cbar = False
        else:
            c = [[0.5, 0.5, 0.5]] * n
            show_cbar = False

    sc = ax.scatter(x_v, y_v, c=c, marker=marker, alpha=alpha,
                    s=10, linewidths=0, zorder=2)
    if show_cbar:
        plt.colorbar(sc, ax=ax)


def _histogram_panel(ax, sio2, y, xlim, ylim, n_bins, bin_scale):
    """Draw 2-D density heatmap."""
    valid = np.isfinite(sio2) & np.isfinite(y) & (y >= 0)
    if not valid.any():
        return

    x_v, y_v = sio2[valid], y[valid]

    xedges = np.linspace(xlim[0], xlim[1], n_bins + 1)
    yedges = np.linspace(ylim[0], ylim[1], n_bins + 1)
    counts, xed, yed = np.histogram2d(x_v, y_v, bins=[xedges, yedges])

    # counts shape is (n_bins, n_bins); imshow expects (rows=y, cols=x)
    img = counts.T

    if bin_scale == 'log':
        with np.errstate(divide='ignore'):
            img = np.log10(np.where(img > 0, img, np.nan))
        vmin, vmax = -0.1, np.nanmax(img)
        cbar_label = 'log$_{10}$(N)'
    else:
        vmin, vmax = 0, np.nanmax(img)
        cbar_label = 'N'

    extent = [xlim[0], xlim[1], ylim[0], ylim[1]]
    im = ax.imshow(
        img, origin='lower', extent=extent, aspect='auto',
        vmin=vmin, vmax=vmax, interpolation='nearest',
    )
    cb = plt.colorbar(im, ax=ax)
    cb.set_label(cbar_label)
