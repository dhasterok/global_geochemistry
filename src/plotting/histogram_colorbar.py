import numpy as np
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable


def histogram_colorbar(mappable, data, cax=None, bins=64,
                       orientation='vertical', size='15%', pad=0.08,
                       clip_out_of_range=False, density=True):
    """Histogram colored by the mappable's cmap+norm, positioned like a colorbar.

    Parameters
    ----------
    mappable : ScalarMappable
        The scatter/image/mesh returned by ax.scatter / imshow / pcolormesh.
    data : array-like
        The values to histogram (typically the same values used for colour).
    cax : Axes, optional
        Axes to draw into.  If None, one is appended next to the mappable's axes
        using ``make_axes_locatable``.
    bins : int
        Number of histogram bins.  Default 64.
    orientation : {'vertical', 'horizontal'}
        Vertical places the histogram to the right of the map (value axis = y).
        Horizontal places it below (value axis = x).
    size : str
        Width (vertical) or height (horizontal) of the appended axes, as a
        percentage string passed to ``make_axes_locatable``.
    pad : float
        Padding between the map axes and the histogram axes.
    clip_out_of_range : bool
        If True, values outside [vmin, vmax] are clipped into the end bins
        rather than discarded.  Default False (out-of-range values are excluded
        so they do not create artificial spikes at the colour limits).
    density : bool
        If True, the histogram is normalised to a probability density.

    Returns
    -------
    cax : matplotlib.axes.Axes
    """
    ax = mappable.axes
    cmap = mappable.get_cmap()
    norm = mappable.norm
    logscale = isinstance(norm, LogNorm)
    vmin, vmax = norm.vmin, norm.vmax

    if cax is None:
        divider = make_axes_locatable(ax)
        side = 'right' if orientation == 'vertical' else 'bottom'
        cax = divider.append_axes(side, size=size, pad=pad)

    x = np.asarray(data, dtype=float)
    x = x[np.isfinite(x)]
    if logscale:
        x = x[x > 0]
        edges = np.geomspace(vmin, vmax, bins + 1)
    else:
        edges = np.linspace(vmin, vmax, bins + 1)

    if clip_out_of_range:
        x = np.clip(x, vmin, vmax)
    counts, _ = np.histogram(x, bins=edges, density=density)

    # Fine colour gradient for the filled rectangle
    ge = np.geomspace(vmin, vmax, 513) if logscale else np.linspace(vmin, vmax, 513)
    gc = np.sqrt(ge[:-1] * ge[1:]) if logscale else 0.5 * (ge[:-1] + ge[1:])

    if orientation == 'vertical':
        # Step-outline polygon (value axis = y)
        ys = np.repeat(edges, 2)[1:-1]      # 2*bins values
        xs = np.repeat(counts, 2)           # 2*bins values
        outline = cax.fill_betweenx(ys, 0, xs, facecolor='none')

        qm = cax.pcolormesh(
            [0, counts.max() * 1.05], ge, gc[:, None],
            cmap=cmap, norm=norm, shading='flat',
        )
        qm.set_clip_path(outline.get_paths()[0], transform=cax.transData)

        if logscale:
            cax.set_yscale('log')
        cax.set_ylim(vmin, vmax)
        # Invert x so the zero baseline is on the right, same side as the
        # y-axis tick labels — bars extend leftward into the histogram axes.
        cax.set_xlim(counts.max() * 1.05, 0)
        cax.yaxis.tick_right()
        cax.yaxis.set_label_position('right')
        cax.set_xticks([])
        for s in ('top', 'bottom', 'left'):
            cax.spines[s].set_visible(False)

    else:
        # Step-outline polygon (value axis = x)
        xs = np.repeat(edges, 2)[1:-1]      # 2*bins values
        ys = np.repeat(counts, 2)           # 2*bins values
        outline = cax.fill_between(xs, 0, ys, facecolor='none')

        qm = cax.pcolormesh(
            ge, [0, counts.max() * 1.05], gc[None, :],
            cmap=cmap, norm=norm, shading='flat',
        )
        qm.set_clip_path(outline.get_paths()[0], transform=cax.transData)

        if logscale:
            cax.set_xscale('log')
        cax.set_xlim(vmin, vmax)
        cax.set_ylim(0, counts.max() * 1.05)
        cax.set_yticks([])
        for s in ('top', 'left', 'right'):
            cax.spines[s].set_visible(False)

    return cax
