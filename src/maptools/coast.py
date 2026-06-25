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
``'geographic'``  (default)
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

from .projections import (cassini_fwd, pyproj_fwd, spilhaus_fwd, waterman_fwd,
                          _ell_int_5)


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

def _geographic(lon, lat, lon0, _lat0):
    # Wrap to [-180, 180] so the result is always in that range regardless of lon0.
    # Without this, lon0=180 gives x ∈ [-360, 0] while the axes expect [-180, 180].
    x = ((np.asarray(lon, dtype=float) - lon0) + 180.0) % 360.0 - 180.0
    return x, np.asarray(lat, dtype=float)


def _make_pyproj(proj_name, lat1=None, lat2=None):
    """Return a forward-projection function for a named PROJ projection.

    For conic projections (lcc, aea) the standard parallels *lat1* and *lat2*
    are forwarded into the PROJ string when provided.
    """
    def _fwd(lon, lat, lon0, lat0, _lat1=lat1, _lat2=lat2):
        parts = [f'+proj={proj_name}', f'+lon_0={lon0}', f'+lat_0={lat0}',
                 '+R=6371000', '+units=m', '+no_defs']
        if _lat1 is not None:
            parts.append(f'+lat_1={_lat1}')
        if _lat2 is not None:
            parts.append(f'+lat_2={_lat2}')
        return pyproj_fwd(lon, lat, ' '.join(parts))
    return _fwd


_PROJECTION_FWD = {
    'geographic':  _geographic,
    'cassini':     lambda lon, lat, lon0, lat0: cassini_fwd(lon, lat, lon0, lat0),
    'spilhaus':    lambda lon, lat, lon0, lat0: spilhaus_fwd(lon, lat),
    'waterman':    lambda lon, lat, lon0, lat0: waterman_fwd(lon, lat, lon0=lon0, style='butterfly'),
    'waterman_m':  lambda lon, lat, lon0, lat0: waterman_fwd(lon, lat, lon0=lon0, style='waterman_m'),
    'mercator':    _make_pyproj('merc'),
    'mollweide':   _make_pyproj('moll'),
    'robinson':    _make_pyproj('robin'),
    'ortho':       _make_pyproj('ortho'),
    'laea':        _make_pyproj('laea'),
    # Conic projections — standard parallels default to ±20° and ±60° for
    # northern (lcc_nh / aea_nh) or southern (lcc_sh / aea_sh) hemisphere.
    # Use 'lcc' / 'aea' directly and pass lat1=/lat2= via plotcoast kwargs.
    'lcc':         _make_pyproj('lcc', lat1=20.0, lat2=60.0),
    'lcc_sh':      _make_pyproj('lcc', lat1=-60.0, lat2=-20.0),
    'aea':         _make_pyproj('aea', lat1=20.0, lat2=60.0),
    'aea_sh':      _make_pyproj('aea', lat1=-60.0, lat2=-20.0),
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
# Pre-projection geographic splitting
# ---------------------------------------------------------------------------

# Octant fold meridians for Waterman projections.  Every polyline segment that
# crosses one of these longitudes must be split there (with a NaN separator)
# before projection so that each sub-segment lies entirely within one octant.
_WATERMAN_FOLDS = np.array([-180.0, -90.0, 0.0, 90.0])


def _presplit_at_folds(
    lon: np.ndarray,
    lat: np.ndarray,
    fold_lons: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Split polyline segments at fold meridians before projection.

    For each finite segment that crosses a meridian in *fold_lons*, inserts the
    linearly interpolated crossing point followed by a NaN separator.  The
    result guarantees every sub-segment lies within a single longitude band.

    Antimeridian crossings (|Δlon| > 180°) are handled separately: the segment
    is split at lon = ±180 and a NaN inserted, preventing wrap-around lines.
    """
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    folds = np.asarray(sorted(fold_lons), dtype=float)
    inner = folds[(folds > -180.0) & (folds < 180.0)]

    out_lon: list[np.ndarray] = [lon[:1]]
    out_lat: list[np.ndarray] = [lat[:1]]

    for i in range(len(lon) - 1):
        lo0, la0 = lon[i], lat[i]
        lo1, la1 = lon[i + 1], lat[i + 1]

        if not (np.isfinite(lo0) and np.isfinite(la0)
                and np.isfinite(lo1) and np.isfinite(la1)):
            out_lon.append(lon[i + 1 : i + 2])
            out_lat.append(lat[i + 1 : i + 2])
            continue

        dlo = lo1 - lo0

        if abs(dlo) > 180.0:
            # Antimeridian crossing — split at ±180
            if dlo < 0:          # going right → left across +180
                t  = (180.0 - lo0) / (lo1 + 360.0 - lo0)
                lc = la0 + t * (la1 - la0)
                out_lon.append(np.array([ 180.0, np.nan]))
                out_lat.append(np.array([   lc,  np.nan]))
            else:                # going left → right across −180
                t  = (-180.0 - lo0) / (lo1 - 360.0 - lo0)
                lc = la0 + t * (la1 - la0)
                out_lon.append(np.array([-180.0, np.nan]))
                out_lat.append(np.array([   lc,  np.nan]))
            out_lon.append(lon[i + 1 : i + 2])
            out_lat.append(lat[i + 1 : i + 2])
            continue

        # Find inner folds between lo0 (exclusive) and lo1 (inclusive).
        # Right-closed so that a segment endpoint landing exactly on a fold
        # meridian is also split (that point belongs to the next octant).
        lo_a, lo_b = min(lo0, lo1), max(lo0, lo1)
        hit = inner[(inner > lo_a) & (inner <= lo_b)]
        if lo1 < lo0:
            hit = hit[::-1]      # traverse in direction lo0 → lo1

        if hit.size == 0:
            out_lon.append(lon[i + 1 : i + 2])
            out_lat.append(lat[i + 1 : i + 2])
            continue

        # For each fold crossing, insert an ε-offset boundary point so that
        # the endpoint stays in the *current* octant, then a NaN separator.
        # Without the ε-offset the fold longitude itself maps to the *next*
        # octant, creating a cross-panel jump in projected space.
        _EPS = 1e-9
        direction = 1 if lo1 >= lo0 else -1   # +1 left→right, −1 right→left
        prev_lo, prev_la = lo0, la0
        for fold_lo in hit:
            t  = (fold_lo - prev_lo) / (lo1 - prev_lo)
            lc = prev_la + t * (la1 - prev_la)
            boundary_lo = fold_lo - direction * _EPS   # stays in current octant
            out_lon.append(np.array([boundary_lo, np.nan]))
            out_lat.append(np.array([lc,          np.nan]))
            prev_lo, prev_la = fold_lo, lc

        out_lon.append(lon[i + 1 : i + 2])
        out_lat.append(lat[i + 1 : i + 2])

    return np.concatenate(out_lon), np.concatenate(out_lat)


def _presplit_at_equator(
    lon: np.ndarray,
    lat: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Split polyline segments that cross the equator (lat = 0).

    The Waterman octahedral projection reflects southern hemisphere points
    across x = altitude in face space (``xMj = 2·alt − xMj``).  For
    longitudes in the *curved* ELD region (|lonr| > _W_lonDiag), the
    equatorial boundary lies at xf < alt, so the north/south faces do not
    meet seamlessly — a structural gap exists in the unfolding.  Splitting
    at lat = 0 replaces the visible cross-equator line with a proper NaN
    separator so the gap reads as a cut, not a spurious connection.

    For longitudes in the *straight* region the boundary is at xf = alt and
    both sides project identically, so the NaN is cosmetically invisible.
    """
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)

    out_lon: list[np.ndarray] = [lon[:1]]
    out_lat: list[np.ndarray] = [lat[:1]]

    for i in range(len(lon) - 1):
        lo0, la0 = lon[i], lat[i]
        lo1, la1 = lon[i + 1], lat[i + 1]

        if not (np.isfinite(lo0) and np.isfinite(la0)
                and np.isfinite(lo1) and np.isfinite(la1)):
            out_lon.append(lon[i + 1 : i + 2])
            out_lat.append(lat[i + 1 : i + 2])
            continue

        # The Waterman projection uses ``south = (lat_rad < 0)`` to decide
        # whether to apply the equatorial reflection.  Any segment that
        # crosses this boundary — i.e. where (la0 < 0) differs from (la1 < 0)
        # — will jump in projected space.  Insert a NaN to separate the two
        # hemisphere segments.  No crossing point is inserted: in the curved
        # ELD region the two octant faces do not meet at a common point, so
        # the gap is intentional (natural cut in the Waterman unfolding).
        if (la0 < 0.0) != (la1 < 0.0):
            out_lon.append(np.array([np.nan]))
            out_lat.append(np.array([np.nan]))

        out_lon.append(lon[i + 1 : i + 2])
        out_lat.append(lat[i + 1 : i + 2])

    return np.concatenate(out_lon), np.concatenate(out_lat)


def _preprocess_lines(
    lon: np.ndarray,
    lat: np.ndarray,
    projection: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Apply projection-specific geographic pre-processing before projection."""
    key = projection.strip().lower()
    if key in ('waterman', 'waterman_m'):
        lon, lat = _presplit_at_folds(lon, lat, _WATERMAN_FOLDS)
        lon, lat = _presplit_at_equator(lon, lat)
        return lon, lat
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


def _clip_geographic(
    x: np.ndarray,
    y: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Clip geographic-projection lines at the ±180 boundary.

    When a polyline segment crosses x = +180 (right edge) or x = −180 (left
    edge), the interpolated crossing point is inserted on both edges so the
    line reaches the frame cleanly rather than stopping short.  NaN separators
    in the input (existing segment boundaries) are preserved unchanged.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 2:
        return x, y

    xs: list[np.ndarray] = [x[:1]]
    ys: list[np.ndarray] = [y[:1]]

    for i in range(len(x) - 1):
        x0, y0 = x[i], y[i]
        x1, y1 = x[i + 1], y[i + 1]

        if not (np.isfinite(x0) and np.isfinite(y0)
                and np.isfinite(x1) and np.isfinite(y1)):
            xs.append(x[i + 1 : i + 2])
            ys.append(y[i + 1 : i + 2])
            continue

        dx = x1 - x0

        if dx < -180.0:
            # Crosses the right edge (+180): exit at +180, re-enter at −180
            t  = (180.0 - x0) / (x1 + 360.0 - x0)
            yc = y0 + t * (y1 - y0)
            xs.append(np.array([ 180.0, np.nan, -180.0, x1]))
            ys.append(np.array([   yc,  np.nan,    yc,  y1]))
        elif dx > 180.0:
            # Crosses the left edge (−180): exit at −180, re-enter at +180
            t  = (-180.0 - x0) / (x1 - 360.0 - x0)
            yc = y0 + t * (y1 - y0)
            xs.append(np.array([-180.0, np.nan, 180.0, x1]))
            ys.append(np.array([   yc,  np.nan,   yc,  y1]))
        else:
            xs.append(x[i + 1 : i + 2])
            ys.append(y[i + 1 : i + 2])

    return np.concatenate(xs), np.concatenate(ys)


def _clip_lines(
    x: np.ndarray,
    y: np.ndarray,
    projection: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Dispatch line clipping to the projection-appropriate function."""
    if projection.strip().lower() == 'geographic':
        return _clip_geographic(x, y)
    return _clip_discontinuities(x, y)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def plotcoast(
    projection: str = 'geographic',
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
        Default ``'geographic'``.
    lon0 : float
        Central meridian in degrees (default 0).  Shifts the map so that
        *lon0* is at the left edge for geographic, or at the centre for
        other projections.
    lat0 : float
        Central parallel in degrees (default 0).  Ignored by geographic.
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
    >>> ax = plotcoast()                              # geographic, global
    >>> ax = plotcoast(lon0=180)                      # Pacific-centred
    >>> ax = plotcoast(projection='cassini', lon0=-150)  # subduction view
    """
    lon, lat = _load_coast()
    lon, lat = _preprocess_lines(lon, lat, projection)
    x, y = _project(lon, lat, projection, lon0, lat0)
    x, y = _clip_lines(x, y, projection)

    if ax is None:
        _, ax = plt.subplots()

    kw.setdefault('color', color)
    kw.setdefault('linewidth', linewidth)
    ax.plot(x, y, **kw)

    ax.set_aspect('equal')
    _set_axis_limits(ax, x, y, projection)
    _draw_mapframe(ax, projection, lon0, lat0,
                   color=color, linewidth=max(linewidth, 0.5))
    _style_map_axes(ax)
    return ax


def geoplot(
    lon,
    lat,
    projection: str = 'geographic',
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
        Projection name (default ``'geographic'``).
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
    lon, lat = _preprocess_lines(lon, lat, projection)
    x, y = _project(lon, lat, projection, lon0, lat0)
    x, y = _clip_lines(x, y, projection)

    if ax is None:
        _, ax = plt.subplots()
    ax.plot(x, y, **kw)
    return ax


def tissot(
    dlat: float = 30.0,
    dlon: float = 30.0,
    radius: float = 5.0,
    n: int = 72,
    projection: str = 'geographic',
    lon0: float = 0.0,
    lat0: float = 0.0,
    ax: plt.Axes | None = None,
    **kw,
) -> plt.Axes:
    """Draw Tissot's indicatrix on the current map.

    Places small circles of angular radius *radius* degrees on a regular
    lat/lon grid.  In an undistorted projection the circles appear circular;
    shape or area distortion deforms them into ellipses, revealing the
    local distortion character of the projection.

    Parameters
    ----------
    dlat, dlon : float
        Grid spacing in degrees (default 30).
    radius : float
        Angular radius of each indicatrix circle in degrees (default 5).
    n : int
        Number of points per circle (default 72).
    projection : str
        Same projection names as :func:`plotcoast`.
    lon0, lat0 : float
        Central meridian and parallel (degrees).
    ax : matplotlib.axes.Axes, optional
        Target axes.  A new figure is created if *None*.
    **kw
        Forwarded to :func:`~matplotlib.axes.Axes.plot`.
        Default: ``color='C0'``, ``linewidth=0.5``, ``alpha=0.6``.

    Returns
    -------
    matplotlib.axes.Axes
    """
    kw.setdefault('color', 'C0')
    kw.setdefault('linewidth', 0.5)
    kw.setdefault('alpha', 0.6)

    if ax is None:
        _, ax = plt.subplots()

    lats = np.arange(-90 + dlat, 90, dlat)
    lons = np.arange(-180 + dlon / 2, 180, dlon)

    r_rad = np.radians(radius)
    theta = np.linspace(0, 2 * np.pi, n, endpoint=False)

    # Build all circles as NaN-separated arrays for a single plot call
    all_lon: list[np.ndarray] = []
    all_lat: list[np.ndarray] = []

    for clat in lats:
        phi0 = np.radians(clat)
        for clon in lons:
            lam0 = np.radians(clon)
            # Parametric small circle via spherical trig
            # Point P at angular distance r_rad from centre (phi0, lam0)
            # bearing theta:
            #   sin(phi) = sin(phi0)*cos(r) + cos(phi0)*sin(r)*cos(az)
            #   lam = lam0 + atan2(sin(az)*sin(r)*cos(phi0),
            #                      cos(r) - sin(phi0)*sin(phi))
            sin_r = np.sin(r_rad)
            cos_r = np.cos(r_rad)
            sin_p0 = np.sin(phi0)
            cos_p0 = np.cos(phi0)

            sin_phi = sin_p0 * cos_r + cos_p0 * sin_r * np.cos(theta)
            sin_phi = np.clip(sin_phi, -1.0, 1.0)
            phi = np.arcsin(sin_phi)
            dlam = np.arctan2(
                np.sin(theta) * sin_r * cos_p0,
                cos_r - sin_p0 * sin_phi,
            )
            lam = lam0 + dlam

            clon_arr = np.degrees(lam)
            clat_arr = np.degrees(phi)

            # Wrap to [-180, 180]
            clon_arr = (clon_arr + 180) % 360 - 180

            # Close the circle
            clon_arr = np.append(clon_arr, clon_arr[0])
            clat_arr = np.append(clat_arr, clat_arr[0])

            all_lon.append(clon_arr)
            all_lon.append(np.array([np.nan]))
            all_lat.append(clat_arr)
            all_lat.append(np.array([np.nan]))

    if not all_lon:
        return ax

    lon_all = np.concatenate(all_lon)
    lat_all = np.concatenate(all_lat)
    lon_all, lat_all = _preprocess_lines(lon_all, lat_all, projection)

    x, y = _project(lon_all, lat_all, projection, lon0, lat0)
    x, y = _clip_lines(x, y, projection)
    ax.plot(x, y, **kw)
    return ax


def geoscatter(
    lon,
    lat,
    projection: str = 'geographic',
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
        Projection name (default ``'geographic'``).
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

def _set_axis_limits(ax, x, y, projection):
    """Set sensible axis limits for the given projection."""
    key = projection.strip().lower()
    if key == 'geographic':
        ax.set_xlim(-180, 180)
        ax.set_ylim(-90, 90)
    elif key == 'cassini':
        from .projections import _R_EARTH
        R = _R_EARTH
        ax.set_xlim(-R * np.pi / 2,  R * np.pi / 2)
        ax.set_ylim(-R * np.pi,      R * np.pi)
    else:
        finite = np.isfinite(x) & np.isfinite(y)
        if finite.any():
            pad_x = 0.02 * (np.nanmax(x) - np.nanmin(x))
            pad_y = 0.02 * (np.nanmax(y) - np.nanmin(y))
            ax.set_xlim(np.nanmin(x) - pad_x, np.nanmax(x) + pad_x)
            ax.set_ylim(np.nanmin(y) - pad_y, np.nanmax(y) + pad_y)


def _frame_path(projection, lon0, lat0, n=360):
    """Return (x, y) arrays that trace the outer boundary of the projection.

    Each projection has a characteristic outline shape:
    - Plate Carrée / Cassini: rectangle (geographic edges projected)
    - pyproj-backed (Mollweide, Robinson, …): the outer oval from lon=±180
    - Spilhaus: analytic diamond from the Adams elliptic-integral scale
    - Waterman: the four octant-boundary meridians forming the butterfly/M outline
    """
    key = projection.strip().lower()
    m = max(n // 4, 90)
    t = np.linspace(0, 1, m, endpoint=False)

    if key == 'geographic':
        # _geographic wraps all x to [-180, 180] regardless of lon0, so the
        # frame is always this fixed rectangle in projected space.
        vx = np.array([-180.0,  180.0,  180.0, -180.0, -180.0])
        vy = np.array([ -90.0,  -90.0,   90.0,   90.0,  -90.0])
        return vx, vy

    elif key == 'spilhaus':
        # Adams World-in-a-Square II diamond boundary.
        # K = complete elliptic integral K(k²=0.5) ≈ 1.854
        # After the 45° rotation in _adams_diamond the diamond vertices
        # lie at (±K√2, 0) and (0, ±K√2).
        K = float(_ell_int_5(np.pi / 2))
        A = K * np.sqrt(2)
        vx = np.array([ A,  0, -A,  0,  A])
        vy = np.array([ 0,  A,  0, -A,  0])
        xs, ys = [], []
        for i in range(4):
            xs.append(vx[i] + (vx[i + 1] - vx[i]) * t)
            ys.append(vy[i] + (vy[i + 1] - vy[i]) * t)
        x = np.append(np.concatenate(xs), vx[0])
        y = np.append(np.concatenate(ys), vy[0])
        return x, y

    elif key == 'cassini':
        # After y-wrapping in cassini_fwd the full strip is always
        # x ∈ [−R·π/2, R·π/2],  y ∈ [−R·π, R·π]  (lat0 shifts the centre
        # latitude but not the frame extents).
        from .projections import _R_EARTH
        R = _R_EARTH
        xc = R * np.pi / 2
        yc = R * np.pi
        vx = np.array([-xc,  xc,  xc, -xc, -xc])
        vy = np.array([-yc, -yc,  yc,  yc, -yc])
        xs, ys = [], []
        for i in range(4):
            xs.append(vx[i] + (vx[i + 1] - vx[i]) * t)
            ys.append(vy[i] + (vy[i + 1] - vy[i]) * t)
        x = np.append(np.concatenate(xs), vx[0])
        y = np.append(np.concatenate(ys), vy[0])
        return x, y

    elif key in ('waterman', 'waterman_m'):
        # Ported from panelborders.m (Waterman MATLAB codebase).
        # Each border segment is a V-path: down one meridian to the south
        # pole then back up the adjacent meridian.  ε-offsets on the fold
        # longitudes ensure OctahedralProjection assigns each point to the
        # correct octant.  The equatorial rows use shifted longitudes
        # (≈ ±18.64°) matching the octahedral face geometry.
        #
        # Latitude breakpoints are the octant boundary latitudes; they must
        # be explicit so the projected path follows the correct octant edges.
        lat_down = np.array([
             89.999999,  71.38865281,  18.7478136,
              0.000001,   0.0,         -0.000001,
            -18.7478136, -71.38865281, -89.999999,
        ])
        lat_up = lat_down[::-1]                     # -89.999999 … +89.999999

        # blon columns from panelborders.m: 4 columns × 9 rows each for the
        # "down" half (top half: rows 1-9 in MATLAB 1-indexing) and the
        # "up" half (bottom half: rows 10-18).
        blon_down = np.array([
            [-179.999999,  -89.999999,   0.000001,  90.000001],
            [-179.999999,  -90.0,         0.0,       90.0     ],
            [-179.999999,  -90.0,         0.0,       90.0     ],
            [-179.999999,  -89.999999,   0.000001,  90.0     ],
            [-161.3603887, -71.36038869, 18.63961131, 108.6396113],
            [-179.999999,  -89.999999,   0.000001,  90.000001],
            [-179.999999,  -89.999999,   0.000001,  90.000001],
            [-179.999999,  -89.999999,   0.000001,  90.000001],
            [-179.999999,  -89.999999,   0.000001,  90.000001],
        ])  # shape (9, 4)

        blon_up = np.array([
            [-90.000001,   -0.000001,  89.999999, 179.999999],
            [-90.000001,   -0.000001,  89.999999, 179.999999],
            [-90.000001,   -0.000001,  89.999999, 179.999999],
            [-90.000001,   -0.000001,  89.999999, 179.999999],
            [-108.6396113, -18.63961131, 71.36038869, 161.3603887],
            [-90.000001,   -0.000001,  89.999999, 179.999999],
            [-90.000001,   -0.000001,  89.999999, 179.999999],
            [-90.000001,   -0.000001,  89.999999, 179.999999],
            [-90.000001,   -0.000001,  89.999999, 179.999999],
        ])  # shape (9, 4)

        xsegs, ysegs = [], []
        for col in range(4):
            v_lon = np.concatenate([blon_down[:, col], blon_up[:, col]])
            v_lat = np.concatenate([lat_down, lat_up])
            xb, yb = _project(v_lon, v_lat, key, lon0, lat0)
            xsegs.extend([xb, np.array([np.nan])])
            ysegs.extend([yb, np.array([np.nan])])
        # No discontinuity clipping: the frame path is deliberately structured
        # and the equatorial kink points would be falsely trimmed at 10%.
        return np.concatenate(xsegs), np.concatenate(ysegs)

    else:
        # Generic: project the geographic boundary rectangle.
        # For pyproj projections (Mollweide, Robinson, …) the lon=±180
        # meridians form the outer oval.
        lons = np.concatenate([
            -180 + 360 * t,              # bottom: left → right
            np.full(m, 180.0),           # right: bottom → top
            180 - 360 * t,               # top: right → left
            np.full(m, -180.0),          # left: top → bottom
        ])
        lats = np.concatenate([
            np.full(m, -90.0),
            -90 + 180 * t,
            np.full(m,  90.0),
            90 - 180 * t,
        ])
        x, y = _project(lons, lats, key, lon0, lat0)
        x, y = _clip_discontinuities(x, y)
        return x, y


_CONIC_KEYS = ('lcc', 'lcc_sh', 'aea', 'aea_sh')


def _draw_mapframe(ax, projection, lon0, lat0, color='k', linewidth=0.5):
    """Draw the projection boundary as a closed frame on *ax*.

    Conic projections (LCC, Albers) do not have a well-defined global frame
    (the poles cannot be projected) so no frame is drawn for them.  Users
    should call ``ax.set_xlim`` / ``ax.set_ylim`` to define the map extent.
    """
    if projection.strip().lower() in _CONIC_KEYS:
        return
    x, y = _frame_path(projection, lon0, lat0)
    ax.plot(x, y, color=color, linewidth=linewidth,
            clip_on=False, zorder=5, solid_capstyle='projecting')


def _style_map_axes(ax):
    """Hide default spines, ticks, and tick labels; keep a data-space frame."""
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
