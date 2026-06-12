"""
Geographic coastline and data plotting.

``plotcoast(projection, lon0, lat0, color, linewidth, ax, **kw)``
    Plot world coastlines on a new or existing axes.  Returns the axes
    handle so the caller can overlay data.

``geoplot(lon, lat, projection, lon0, lat0, ax, **kw)``
    Project geographic line data and plot with automatic discontinuity
    clipping (no cross-map artefacts from dateline/boundary crossings).

``geoscatter(lon, lat, projection, lon0, lat0, ax, **kw)``
    Project and scatter-plot geographic point data.

Supported projections
---------------------
``'platecarree'``  (default)
    Simple equirectangular — x = lon − lon0, y = lat.
``'cassini'``
    Spherical Cassini-Soldner.  Maps the full globe onto a vertical strip;
    setting lon0 to the Pacific (≈ −150°) produces a "south-pole to
    south-pole" view that is excellent for global subduction.
``'mercator'``, ``'mollweide'``, ``'robinson'``, ``'ortho'``, ``'laea'``
    Standard projections via pyproj (already in environment).  lon0 is
    forwarded as +lon_0 and lat0 as +lat_0.

Phase-B additions (not yet implemented)
----------------------------------------
``'spilhaus'``   — ocean-centric (Adams/Spilhaus); port from ``adams.m``.
``'waterman'``   — Waterman butterfly; port from ``OctahedralProjection.m``.
``'waterman_m'`` — Waterman M style.

Coast data
----------
Loaded from ``src/data/coast.npz``.  If that file is missing the module
auto-converts ``matlab/maptools/coast.mat`` (via ``scipy.io.loadmat``) and
saves the NPZ for subsequent calls.
"""

from __future__ import annotations

import pathlib

import numpy as np
import matplotlib.pyplot as plt

from .projections import cassini_fwd, pyproj_fwd, spilhaus_fwd, waterman_fwd


# ---------------------------------------------------------------------------
# Data paths
# ---------------------------------------------------------------------------

_HERE     = pathlib.Path(__file__).resolve().parent
_SRC_DATA = _HERE.parent / 'data'
_COAST_NPZ = _SRC_DATA / 'coast.npz'
_COAST_MAT = _HERE.parent.parent / 'matlab' / 'maptools' / 'coast.mat'


# ---------------------------------------------------------------------------
# Projection registry
# ---------------------------------------------------------------------------

def _platecarree(lon, lat, lon0, lat0):
    return (np.asarray(lon, dtype=float) - lon0,
            np.asarray(lat, dtype=float))


def _make_pyproj(proj_name):
    """Return a forward-projection function for a named PROJ projection."""
    def _fwd(lon, lat, lon0, lat0):
        proj_str = (f'+proj={proj_name} +lon_0={lon0} +lat_0={lat0}'
                    f' +R=6371000 +units=m +no_defs')
        return pyproj_fwd(lon, lat, proj_str)
    return _fwd


_PROJECTION_FWD = {
    'platecarree': _platecarree,
    'cassini':     lambda lon, lat, lon0, lat0: cassini_fwd(lon, lat, lon0, lat0),
    'spilhaus':    lambda lon, lat, lon0, lat0: spilhaus_fwd(lon, lat),
    'waterman':    lambda lon, lat, lon0, lat0: waterman_fwd(lon, lat, lon0=lon0, style='butterfly'),
    'waterman_m':  lambda lon, lat, lon0, lat0: waterman_fwd(lon, lat, lon0=lon0, style='waterman_m'),
    'mercator':    _make_pyproj('merc'),
    'mollweide':   _make_pyproj('moll'),
    'robinson':    _make_pyproj('robin'),
    'ortho':       _make_pyproj('ortho'),
    'laea':        _make_pyproj('laea'),
}


# ---------------------------------------------------------------------------
# Coast data loading
# ---------------------------------------------------------------------------

def _load_coast() -> tuple[np.ndarray, np.ndarray]:
    """Return (lon, lat) coastline arrays with NaN segment separators."""
    if _COAST_NPZ.exists():
        d = np.load(_COAST_NPZ)
        return d['lon'], d['lat']

    # Auto-convert from MATLAB .mat on first use
    if not _COAST_MAT.exists():
        raise FileNotFoundError(
            f'Coast data not found.  Expected either:\n'
            f'  {_COAST_NPZ}\n'
            f'  {_COAST_MAT}'
        )
    from scipy.io import loadmat
    d = loadmat(str(_COAST_MAT))
    lon = d['clon'].ravel().astype(float)
    lat = d['clat'].ravel().astype(float)
    _COAST_NPZ.parent.mkdir(parents=True, exist_ok=True)
    np.savez(_COAST_NPZ, lon=lon, lat=lat)
    return lon, lat


# ---------------------------------------------------------------------------
# Projection + discontinuity clipping
# ---------------------------------------------------------------------------

def _project(
    lon: np.ndarray,
    lat: np.ndarray,
    projection: str,
    lon0: float,
    lat0: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Apply the named forward projection."""
    key = projection.strip().lower()
    if key not in _PROJECTION_FWD:
        raise ValueError(
            f"Unknown projection '{projection}'.  "
            f"Available: {sorted(_PROJECTION_FWD)}"
        )
    return _PROJECTION_FWD[key](lon, lat, lon0, lat0)


def _clip_discontinuities(
    x: np.ndarray,
    y: np.ndarray,
    threshold_frac: float = 0.1,
) -> tuple[np.ndarray, np.ndarray]:
    """Insert NaN where projected coordinates jump discontinuously.

    Detects edges in a NaN-separated polyline where consecutive finite
    points jump by more than *threshold_frac* × map extent.  These jumps
    are artefacts of features that cross the map boundary (e.g. a coastline
    crossing the antimeridian).

    Parameters
    ----------
    x, y : numpy.ndarray
        Projected coordinate arrays (already contain NaN segment separators).
    threshold_frac : float
        Fraction of total map extent used as the jump threshold (default 0.1).

    Returns
    -------
    x, y : numpy.ndarray
        Arrays with additional NaN values inserted at discontinuities.
    """
    x = x.copy()
    y = y.copy()

    finite = np.isfinite(x) & np.isfinite(y)
    if finite.sum() < 2:
        return x, y

    x_range = np.nanmax(x) - np.nanmin(x)
    y_range = np.nanmax(y) - np.nanmin(y)
    if x_range == 0 and y_range == 0:
        return x, y

    thresh_x = threshold_frac * x_range if x_range > 0 else np.inf
    thresh_y = threshold_frac * y_range if y_range > 0 else np.inf

    dx = np.abs(np.diff(x))
    dy = np.abs(np.diff(y))

    # A jump is a large consecutive difference that is NOT a NaN→finite or
    # finite→NaN transition (those are the existing segment separators).
    both_finite = finite[:-1] & finite[1:]
    jumps = both_finite & ((dx > thresh_x) | (dy > thresh_y))

    # Insert NaN *after* the jump start (i.e. before index i+1)
    insert_pos = np.where(jumps)[0] + 1
    if insert_pos.size:
        x = np.insert(x, insert_pos, np.nan)
        y = np.insert(y, insert_pos, np.nan)
    return x, y


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def plotcoast(
    projection: str = 'platecarree',
    lon0: float = 0.0,
    lat0: float = 0.0,
    color='k',
    linewidth: float = 0.24,
    ax: plt.Axes | None = None,
    **kw,
) -> plt.Axes:
    """Plot world coastlines using the specified projection.

    Parameters
    ----------
    projection : str
        Projection name.  See module docstring for supported values.
        Default ``'platecarree'``.
    lon0 : float
        Central meridian in degrees (default 0).  Shifts the map so that
        *lon0* is at the left edge for plate carrée, or at the centre for
        other projections.
    lat0 : float
        Central parallel in degrees (default 0).  Ignored by plate carrée.
    color : matplotlib color
        Coastline colour (default black ``'k'``).
    linewidth : float
        Line width in points (default 0.24, matching MATLAB original).
    ax : matplotlib.axes.Axes, optional
        Target axes.  A new figure/axes is created if *None*.
    **kw
        Additional keyword arguments forwarded to :func:`matplotlib.axes.Axes.plot`.

    Returns
    -------
    matplotlib.axes.Axes
        The axes used for plotting (new or passed-in).

    Examples
    --------
    >>> ax = plotcoast()                             # plate carrée, global
    >>> ax = plotcoast(lon0=180)                     # Pacific-centred
    >>> ax = plotcoast(projection='cassini', lon0=-150)  # subduction view
    """
    lon, lat = _load_coast()
    x, y = _project(lon, lat, projection, lon0, lat0)
    x, y = _clip_discontinuities(x, y)

    if ax is None:
        _, ax = plt.subplots()

    kw.setdefault('color', color)
    kw.setdefault('linewidth', linewidth)
    ax.plot(x, y, **kw)

    ax.set_aspect('equal')
    _set_axis_limits(ax, x, y, projection, lon0)
    return ax


def geoplot(
    lon,
    lat,
    projection: str = 'platecarree',
    lon0: float = 0.0,
    lat0: float = 0.0,
    ax: plt.Axes | None = None,
    **kw,
) -> plt.Axes:
    """Project and plot geographic line data.

    Projects *lon*/*lat* to map coordinates, clips discontinuities caused by
    features crossing the map boundary, then calls ``ax.plot()``.

    Parameters
    ----------
    lon, lat : array-like
        Geographic coordinates in decimal degrees.  May include NaN segment
        separators if the data is already broken into polylines.
    projection : str
        Projection name (default ``'platecarree'``).
    lon0, lat0 : float
        Central meridian and parallel (degrees).
    ax : matplotlib.axes.Axes, optional
        Target axes.  A new figure is created if *None*.
    **kw
        Forwarded to :func:`~matplotlib.axes.Axes.plot`.

    Returns
    -------
    matplotlib.axes.Axes
    """
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    x, y = _project(lon, lat, projection, lon0, lat0)
    x, y = _clip_discontinuities(x, y)

    if ax is None:
        _, ax = plt.subplots()
    ax.plot(x, y, **kw)
    return ax


def geoscatter(
    lon,
    lat,
    projection: str = 'platecarree',
    lon0: float = 0.0,
    lat0: float = 0.0,
    ax: plt.Axes | None = None,
    **kw,
) -> plt.Axes:
    """Project and scatter-plot geographic point data.

    Parameters
    ----------
    lon, lat : array-like
        Geographic coordinates in decimal degrees.
    projection : str
        Projection name (default ``'platecarree'``).
    lon0, lat0 : float
        Central meridian and parallel (degrees).
    ax : matplotlib.axes.Axes, optional
        Target axes.  A new figure is created if *None*.
    **kw
        Forwarded to :func:`~matplotlib.axes.Axes.scatter`.

    Returns
    -------
    matplotlib.axes.Axes
    """
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    x, y = _project(lon, lat, projection, lon0, lat0)

    if ax is None:
        _, ax = plt.subplots()
    ax.scatter(x, y, **kw)
    return ax


# ---------------------------------------------------------------------------
# Axis helpers
# ---------------------------------------------------------------------------

def _set_axis_limits(ax, x, y, projection, lon0):
    """Set sensible axis limits for the given projection."""
    key = projection.strip().lower()
    if key == 'platecarree':
        ax.set_xlim(lon0 - 180, lon0 + 180)
        ax.set_ylim(-90, 90)
    else:
        finite = np.isfinite(x) & np.isfinite(y)
        if finite.any():
            pad_x = 0.02 * (np.nanmax(x) - np.nanmin(x))
            pad_y = 0.02 * (np.nanmax(y) - np.nanmin(y))
            ax.set_xlim(np.nanmin(x) - pad_x, np.nanmax(x) + pad_x)
            ax.set_ylim(np.nanmin(y) - pad_y, np.nanmax(y) + pad_y)
