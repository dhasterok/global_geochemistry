"""
Custom and pyproj-backed geographic projections.

All forward functions accept arrays of (lon, lat) in decimal degrees and
return (x, y).  Cassini returns metres; Spilhaus and Waterman return
dimensionless units (the elliptic-integral and octahedral scales
respectively).  Inverse functions are the reverse.

Available
---------
``cassini_fwd / cassini_inv``
    Spherical Cassini-Soldner (R = 6 371 km).  Port of ``cassini.m`` /
    ``cassini_inv.m``.

``spilhaus_fwd``
    Spilhaus ocean-centric projection.  Rotates the globe so pole-points
    land in continental interiors (≈115°E 30°N and 56°W 30°S), then
    applies Adams World-in-a-Square II (diamond conformal).  Port of
    ``adams.m`` (GDrive maptools/).

``waterman_fwd``
    Waterman Butterfly or M projection.  Unfolds the sphere onto a regular
    octahedron and lays the faces flat.  Port of ``OctahedralProjection.m``
    (GDrive waterman/).

``pyproj_fwd / pyproj_inv``
    Thin wrappers around :mod:`pyproj` for any standard PROJ projection.
"""

from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# Cassini-Soldner (spherical)
# ---------------------------------------------------------------------------

_R_EARTH = 6_371_000.0   # metres


def cassini_fwd(
    lon: np.ndarray,
    lat: np.ndarray,
    lon0: float = 0.0,
    lat0: float = 0.0,
    R: float = _R_EARTH,
) -> tuple[np.ndarray, np.ndarray]:
    """Spherical Cassini-Soldner forward projection.

    Maps the full globe onto a vertical strip
    x ∈ [−Rπ/2, Rπ/2], y ∈ [−2R, 2R].  Points in the far hemisphere
    (|lon − lon0| > 90°) wrap around, placing both poles at the top and
    bottom of the strip — useful for "south-pole to south-pole" views of
    global subduction when centred on the Pacific.

    Parameters
    ----------
    lon, lat : array-like
        Geographic coordinates (decimal degrees).  NaN values are preserved.
    lon0 : float
        Central meridian (degrees).  Default 0.
    lat0 : float
        Central parallel (degrees).  Default 0.
    R : float
        Sphere radius in metres.  Default 6 371 000 m.

    Returns
    -------
    x, y : numpy.ndarray
        Projected coordinates in metres.
    """
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    phi  = np.radians(lat)
    lam  = np.radians(lon - lon0)
    phi0 = np.radians(lat0)

    b = np.cos(phi) * np.sin(lam)
    x = R * np.arcsin(np.clip(b, -1.0, 1.0))
    y = R * (np.arctan2(np.tan(phi), np.cos(lam)) - phi0)
    return x, y


def cassini_inv(
    x: np.ndarray,
    y: np.ndarray,
    lon0: float = 0.0,
    lat0: float = 0.0,
    R: float = _R_EARTH,
) -> tuple[np.ndarray, np.ndarray]:
    """Spherical Cassini-Soldner inverse projection.

    Parameters
    ----------
    x, y : array-like
        Projected coordinates in metres.
    lon0, lat0 : float
        Central meridian and parallel (degrees).
    R : float
        Sphere radius in metres.

    Returns
    -------
    lon, lat : numpy.ndarray
        Geographic coordinates in decimal degrees.  Longitude normalised to
        [−180, 180].
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    phi0 = np.radians(lat0)

    xr = x / R
    yr = y / R + phi0

    lat = np.degrees(np.arcsin(np.sin(yr) * np.cos(xr)))
    lon = (np.degrees(np.arctan2(np.tan(xr), np.cos(yr))) + lon0 + 180.0) % 360.0 - 180.0
    return lon, lat


# ---------------------------------------------------------------------------
# Spilhaus ocean-centric projection
# (Adams World-in-a-Square II with globe rotation)
# ---------------------------------------------------------------------------

_ELL_C0 = 2.19174570831038
_ELL_C  = np.array([-8.58691003636495e-07,  2.02692115653689e-07,
                      3.12960480765314e-05,  5.30394739921063e-05,
                     -1.28046446806130e-03, -5.75574836830288e-03,
                      9.14203033408211e-02])


def _ell_int_5(phi: np.ndarray) -> np.ndarray:
    """Elliptic integral of the first kind, k²=0.5 (Chebyshev approximation).

    Precision better than 1e-7.  Port of the MATLAB ``ell_int_5`` sub-function
    in ``adams.m``.
    """
    y  = 2.0 * (phi * (2.0 / np.pi)) ** 2 - 1.0
    y2 = 2.0 * y
    d1 = np.zeros_like(phi)
    d2 = np.zeros_like(phi)
    for ci in _ELL_C:
        d1, d2 = y2 * d1 - d2 + ci, d1
    return phi * (y * d1 - d2 + 0.5 * _ELL_C0)


def _adams_diamond(lon: np.ndarray, lat: np.ndarray,
                   lon0: float = 0.0, lat0: float = 0.0):
    """Adams World-in-a-Square II (diamond / world2) forward projection.

    Vectorised port of the ``"world2"``/``"diamond"`` branch in ``adams.m``.
    Returns dimensionless (x, y); NaN in → NaN out.
    """
    lon = np.asarray(lon, dtype=float).copy()
    lat = np.asarray(lat, dtype=float).copy()
    nan_mask = np.isnan(lon) | np.isnan(lat)

    lon -= lon0
    lon  = (lon + 180.0) % 360.0 - 180.0
    lam  = np.radians(lon)
    phi  = np.radians(lat - lat0)

    spp   = np.tan(0.5 * phi)
    a_raw = np.cos(np.arcsin(np.clip(spp, -1.0, 1.0))) * np.sin(0.5 * lam)

    indm = (spp + a_raw) < 0
    indn = (spp - a_raw) < 0

    b = np.arccos(np.clip(spp,   -1.0, 1.0))
    a = np.arccos(np.clip(a_raw, -1.0, 1.0))

    m = np.arcsin(np.sqrt(np.maximum(0.0, 1.0 + np.minimum(0.0, np.cos(a + b)))))
    n = np.arcsin(np.sqrt(np.abs(  1.0 - np.maximum(0.0, np.cos(a - b)))))
    m[indm] = -m[indm]
    n[indn] = -n[indn]

    x = _ell_int_5(m)
    y = _ell_int_5(n)

    # 45° rotation for diamond layout
    x, y = (x - y) / np.sqrt(2), (x + y) / np.sqrt(2)

    x[nan_mask] = np.nan
    y[nan_mask] = np.nan
    return x, y


def spilhaus_fwd(
    lon: np.ndarray,
    lat: np.ndarray,
    lon0: float = 0.0,   # unused — Spilhaus has a fixed projection centre
    lat0: float = 0.0,   # unused
) -> tuple[np.ndarray, np.ndarray]:
    """Spilhaus ocean-centric forward projection.

    Rotates the globe so the two land-pole points lie at ≈ 115°E 30°N and
    56°W 30°S (leaving the world ocean as a connected body), then applies
    the Adams World-in-a-Square II (diamond) conformal projection.

    Port of the ``'Spilhaus'`` branch in ``adams.m`` (GDrive maptools/).
    Output is dimensionless (roughly ±2.6 in each axis).

    Parameters
    ----------
    lon, lat : array-like
        Geographic coordinates (decimal degrees).
    lon0, lat0 : float
        Ignored — the projection centre is fixed by the globe rotation.

    Returns
    -------
    x, y : numpy.ndarray
        Dimensionless projected coordinates.
    """
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    shape = lon.shape
    lon_r = lon.ravel()
    lat_r = lat.ravel()

    # Globe rotation: Ry(-60°) · Rz(-115°)
    tz = np.radians(-115.0)
    ty = np.radians(-60.0)

    lam = np.radians(lon_r)
    phi = np.radians(lat_r)
    xc  = np.cos(phi) * np.cos(lam)
    yc  = np.cos(phi) * np.sin(lam)
    zc  = np.sin(phi)

    cz, sz = np.cos(tz), np.sin(tz)
    cy, sy = np.cos(ty), np.sin(ty)

    # Apply Rz then Ry
    xr =  cz * xc - sz * yc
    yr =  sz * xc + cz * yc
    zr =  zc
    x2 =  cy * xr + sy * zr
    y2 =  yr
    z2 = -sy * xr + cy * zr

    lon2 = np.degrees(np.arctan2(y2, x2)).reshape(shape)
    lat2 = np.degrees(np.arctan2(z2, np.sqrt(x2**2 + y2**2))).reshape(shape)

    # Adams diamond with internal centre at lon0=-24, lat0=0
    return _adams_diamond(lon2, lat2, lon0=-24.0, lat0=0.0)


# ---------------------------------------------------------------------------
# Waterman Butterfly / M projections
# ---------------------------------------------------------------------------

_S3 = np.sqrt(3)

_WATERMAN_CONFIGS: dict[str, dict] = {
    'butterfly': {
        'altitude':   2 * _S3,
        'fullHeight': 4 / _S3,
        'octants': np.array([
            [0.0, -1/_S3, -np.pi/2, -3*np.pi/4],
            [0.0, -1/_S3, -np.pi/6,   -np.pi/4],
            [0.0, -1/_S3,  np.pi/6,    np.pi/4],
            [0.0, -1/_S3,  np.pi/2,  3*np.pi/4],
        ]),
    },
    'waterman_m': {
        'altitude':   2 * _S3,
        'fullHeight': _S3,
        'octants': np.array([
            [-1.0, 0.0, -np.pi/6, -3*np.pi/4],
            [-1.0, 0.0,  np.pi/6,   -np.pi/4],
            [ 1.0, 0.0, -np.pi/6,    np.pi/4],
            [ 1.0, 0.0,  np.pi/6,  3*np.pi/4],
        ]),
    },
}

# Waterman face geometry constants
_W_dYdL   = np.array([0.0, 2.0/np.pi, 6.0/np.pi, (4.0 + np.sqrt(8.0)) / np.pi])
_W_lonDiag = np.pi / (4.0 + np.sqrt(8.0))
_W_sin15   = (_S3 - 1.0) / np.sqrt(8.0)
_W_cos15   = (_S3 + 1.0) / np.sqrt(8.0)


def _linInterp(x, x1, x2, y1, y2):
    """Linear interpolation; safe against x1==x2."""
    denom = x2 - x1
    t = np.where(np.abs(denom) > 0, (x - x1) / np.where(denom != 0, denom, 1.0), 0.0)
    return y1 + (y2 - y1) * t


def _joint_heights(lon: np.ndarray):
    """Compute the four Equal-Line-Delineation (ELD) keypoints for each meridian.

    Parameters
    ----------
    lon : (N,) array
        Absolute longitude values in radians (0 ≤ lon ≤ π/2).

    Returns
    -------
    xELD, yELD : (N, 4) arrays
        x and y coordinates of the four ELD boundary points.
    """
    N = len(lon)
    xELD = np.tile([(_S3 - 1) / 2, _S3 / 2, 3 * _S3 / 2, 0.0], (N, 1)).astype(float)
    yELD = np.outer(lon, _W_dYdL)

    straight = lon <= _W_lonDiag
    curved   = ~straight

    xELD[straight, 3] = 2.0 * _S3
    # yELD[straight, 3] already = dYdL[3] * lon from outer product

    if curved.any():
        d = lon[curved] - _W_lonDiag
        xELD[curved, 3] = -d * _W_dYdL[3] * _W_sin15 + 2.0 * _S3
        yELD[curved, 3] =  d * _W_dYdL[3] * _W_cos15 + 1.0

    return xELD, yELD


def _face_forward(lat: np.ndarray, lon: np.ndarray):
    """Project a point on a generic octant face.

    Parameters
    ----------
    lat, lon : (N,) arrays
        Absolute latitude and longitude in radians (both ≥ 0).

    Returns
    -------
    x, y : (N,) arrays
        Coordinates in octant space.
    """
    N = len(lat)
    xELD, yELD = _joint_heights(lon)       # (N, 4)

    dx = np.diff(xELD, axis=1)             # (N, 3)
    dy = np.diff(yELD, axis=1)
    seg = np.sqrt(dx**2 + dy**2)           # (N, 3) segment lengths

    totlen = np.sum(seg, axis=1)           # (N,)

    # Distance from pole
    d = _linInterp(lat, np.pi / 2, 0.0, 0.0, totlen)  # (N,)

    # Remaining budget at each segment start: (N, 3)
    cumlen = np.hstack([np.zeros((N, 1)),
                        np.cumsum(seg[:, :2], axis=1)])
    desir = d[:, None] - cumlen

    # Which segment handles the point?  (port of MATLAB i = 4 - sum(ind))
    fits  = seg >= desir                   # (N, 3)
    sidx  = np.clip(3 - np.sum(fits, axis=1), 0, 2)   # 0-based

    row = np.arange(N)
    dl  = desir[row, sidx]
    li  = seg[row, sidx]
    x0  = xELD[row, sidx];     x1 = xELD[row, sidx + 1]
    y0  = yELD[row, sidx];     y1 = yELD[row, sidx + 1]

    x = _linInterp(dl, 0.0, li, x0, x1)
    y = _linInterp(dl, 0.0, li, y0, y1)
    return x, y


def _waterman_core(lat_rad: np.ndarray, lon_rad: np.ndarray, cfg: dict):
    """Core Waterman forward projection (radians → dimensionless x, y)."""
    alt  = cfg['altitude']
    octs = cfg['octants']

    # Assign each point to an octant quadrant (0-based)
    q = np.floor(2.0 * (lon_rad + np.pi) / np.pi).astype(int) % 4
    q = np.clip(q, 0, 3)

    # Longitude relative to the octant's reference meridian
    lonr = (lon_rad - octs[q, 3] + np.pi) % (2.0 * np.pi) - np.pi

    # Octant vertex position and rotation angle
    xP    = octs[q, 0] * alt
    yP    = octs[q, 1] * alt
    theta = octs[q, 2]

    # Project onto the canonical face (lat ≥ 0, lonr ≥ 0)
    xf, yf = _face_forward(np.abs(lat_rad), np.abs(lonr))

    xMj = xf.copy()
    yMj = yf.copy()

    # Reflect southern hemisphere over the equatorial edge
    south = lat_rad < 0
    xMj[south] = 2.0 * alt - xMj[south]

    # Reflect west longitudes
    west = lonr < 0
    yMj[west] = -yMj[west]

    # Rotate onto the octant and translate
    x = xP + np.sin(theta) * xMj + np.cos(theta) * yMj
    y = yP - np.cos(theta) * xMj + np.sin(theta) * yMj \
        + cfg['fullHeight'] * alt / 2.0

    return x, y


def waterman_fwd(
    lon: np.ndarray,
    lat: np.ndarray,
    lon0: float = 0.0,
    lat0: float = 0.0,   # unused — octahedral projections do not support lat0
    style: str = 'butterfly',
) -> tuple[np.ndarray, np.ndarray]:
    """Waterman Butterfly or M projection forward transform.

    Unfolds the sphere onto a regular octahedron and arranges the faces in
    either the classic four-lobe Butterfly layout or the side-by-side M
    layout.  Port of ``OctahedralProjection.m`` (GDrive waterman/).

    Parameters
    ----------
    lon, lat : array-like
        Geographic coordinates (decimal degrees).
    lon0 : float
        Longitudinal shift in degrees (shifts which meridian is centred).
    lat0 : float
        Ignored.
    style : {'butterfly', 'waterman_m'}
        Map layout.

    Returns
    -------
    x, y : numpy.ndarray
        Dimensionless projected coordinates.
    """
    key = style.lower()
    if key not in _WATERMAN_CONFIGS:
        raise ValueError(f"Unknown Waterman style '{style}'. "
                         f"Choose from {sorted(_WATERMAN_CONFIGS)}.")
    cfg = _WATERMAN_CONFIGS[key]

    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    shape = lon.shape
    lon_r = lon.ravel().copy()
    lat_r = lat.ravel()

    if lon0 != 0.0:
        lon_r = (lon_r + lon0 + 180.0) % 360.0 - 180.0

    lon_rad = np.radians(lon_r)
    lat_rad = np.radians(lat_r)

    x = np.full(lon_r.shape, np.nan)
    y = np.full(lon_r.shape, np.nan)

    valid = np.isfinite(lon_rad) & np.isfinite(lat_rad)
    if valid.any():
        x[valid], y[valid] = _waterman_core(lat_rad[valid], lon_rad[valid], cfg)

    return x.reshape(shape), y.reshape(shape)


# ---------------------------------------------------------------------------
# pyproj wrapper for standard PROJ projections
# ---------------------------------------------------------------------------

_PROJ_SENTINEL = 1e20   # pyproj returns ~1e30 for out-of-domain points


def pyproj_fwd(
    lon: np.ndarray,
    lat: np.ndarray,
    proj_string: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Forward projection via pyproj.

    Parameters
    ----------
    lon, lat : array-like
        Geographic coordinates (decimal degrees).
    proj_string : str
        Any PROJ4/PROJ6 projection string, e.g.
        ``'+proj=moll +lon_0=-150 +R=6371000'``.

    Returns
    -------
    x, y : numpy.ndarray
        Projected coordinates in metres.  Out-of-domain points and PROJ
        sentinel values (|v| > 1e20) are replaced with NaN.
    """
    from pyproj import Proj
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    p = Proj(proj_string)
    x, y = p(lon, lat)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x[np.abs(x) > _PROJ_SENTINEL] = np.nan
    y[np.abs(y) > _PROJ_SENTINEL] = np.nan
    return x, y


def pyproj_inv(
    x: np.ndarray,
    y: np.ndarray,
    proj_string: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Inverse projection via pyproj.

    Parameters
    ----------
    x, y : array-like
        Projected coordinates in metres.
    proj_string : str
        PROJ projection string (same as used for the forward transform).

    Returns
    -------
    lon, lat : numpy.ndarray
        Geographic coordinates in decimal degrees.
    """
    from pyproj import Proj
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    p = Proj(proj_string)
    lon, lat = p(x, y, inverse=True)
    return np.asarray(lon, dtype=float), np.asarray(lat, dtype=float)
