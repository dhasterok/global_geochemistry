"""
Geodetic coordinate transforms and great-circle paths.

Coordinate transforms
---------------------
``ll2utm(lon, lat, datum='wgs84')``
    Geographic (lon, lat) → UTM (easting, northing, zone).

``utm2ll(easting, northing, zone, datum='wgs84')``
    UTM (easting, northing, zone) → geographic (lon, lat).

``mercator(lon, lat, datum='wgs84', lon0=0.0)``
    Forward ellipsoidal Mercator projection.

``inv_mercator(x, y, datum='wgs84', lon0=0.0)``
    Inverse Mercator projection.

``polarstereo_fwd(lon, lat, lat_c=-70, lon0=0, a=..., e=...)``
    Polar stereographic forward projection (NSIDC/SCAR convention).

``polarstereo_inv(x, y, lat_c=-70, lon0=0, a=..., e=...)``
    Polar stereographic inverse projection.

``authalic_lat(lat, datum='wgs84')``
    Geodetic → authalic latitude and equivalent sphere radius.

Spherical geometry
------------------
``gcpoints(P, Q, n)``
    Compute *n* equally-spaced points along the great-circle arc from P
    to Q (each a ``(lon, lat)`` pair in degrees).

``sphangle(P, Q)``
    Angular separation between two points on the sphere (degrees).

``smallcircle(lon0, lat0, delta, n=360)``
    Points on the small circle of angular radius *delta* about (lon0, lat0).

Area calculations
-----------------
``pixarea(lat1, lat2, lon1, lon2, datum=None)``
    Fractional area of a lat/lon quad pixel on a sphere.

``rectarea(lonlim, latlim, R=6371.0)``
    Area of a lat/lon rectangle in km².

All longitude/latitude in decimal degrees.  Negative zone values from
:func:`ll2utm` / :func:`utm2ll` indicate the southern hemisphere.

References
----------
Snyder, J.P. (1987) *Map Projections — A Working Manual*.
USGS Professional Paper 1395. US Government Printing Office.
"""

from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# Reference ellipsoid parameters
# ---------------------------------------------------------------------------

_DATUMS: dict[str, tuple[float, float]] = {
    # (semimajor a [m], semiminor b [m])
    'nad27':         (6_378_206.4,   6_356_583.8),
    'clarke1866':    (6_378_206.4,   6_356_583.8),
    'clarke1880':    (6_378_249.1,   6_356_514.9),
    'nad83':         (6_378_137.0,   6_356_752.3),
    'wgs84':         (6_378_137.0,   6_356_752.3),
    'grs80':         (6_378_137.0,   6_356_752.3),
    'wgrs80':        (6_378_137.0,   6_356_752.3),
    'wgrs84':        (6_378_137.0,   6_356_752.3),
    'nztm':          (6_378_137.0,   6_356_752.3),
    'australian1965':(6_378_160.0,   6_356_774.7),
    'australian1984':(6_378_160.0,   6_356_774.7),
    'gda84':         (6_378_160.0,   6_356_774.7),
    'australian1994':(6_377_137.0,   6_355_755.7),
    'gda94':         (6_377_137.0,   6_355_755.7),
    'krasovsky1940': (6_378_245.0,   6_356_863.0),
    'intl1924':      (6_378_388.0,   6_356_911.9),
    'hayford1909':   (6_378_388.0,   6_356_911.9),
    'nzmg':          (6_378_388.0,   6_356_911.9),
    'airy1830':      (6_377_563.4,   6_356_356.9),
    'bessel1841':    (6_377_397.155, 6_356_079.0),
    'everest1830':   (6_377_376.3,   6_356_075.4),
}


def _get_datum(name: str) -> tuple[float, float]:
    key = name.lower()
    if key not in _DATUMS:
        raise ValueError(f"Unknown datum '{name}'. Available: {sorted(_DATUMS)}")
    return _DATUMS[key]


def _eccentricity(a: float, b: float) -> float:
    return np.sqrt(1 - (b / a)**2)


# ---------------------------------------------------------------------------
# UTM false origins by datum
# ---------------------------------------------------------------------------

def _false_origin(datum: str, zone: float | None = None) -> tuple[float, float]:
    """Return (false_easting, false_northing) for a datum."""
    d = datum.lower()
    if d in ('australian1984', 'gda84'):
        return 0.0, 1e7
    if d in ('australian1994', 'gda94'):
        return 5e5, 1e7
    if d == 'nztm':
        return 1_600_000.0, 10_000_000.0
    return 5e5, 0.0


# ---------------------------------------------------------------------------
# ll2utm
# ---------------------------------------------------------------------------

def ll2utm(
    lon: float | np.ndarray,
    lat: float | np.ndarray,
    datum: str = 'wgs84',
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert geographic coordinates to UTM.

    Parameters
    ----------
    lon, lat : array-like
        Longitude and latitude in decimal degrees.
    datum : str
        Reference ellipsoid name (default ``'wgs84'``).

    Returns
    -------
    easting, northing : ndarray
        UTM coordinates in metres.
    zone : ndarray of int
        UTM zone number (1–60).  Negative for southern hemisphere.
    """
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)

    # Wrap longitudes > 180
    lon = np.where(lon >= 180, lon - 360, lon)

    a, b = _get_datum(datum)
    ecc  = _eccentricity(a, b)
    feast, fnorth = _false_origin(datum)

    zone  = (np.floor((lon + 180) / 6) + 1).astype(int)
    lon0  = np.radians(6 * (zone - 1) + 3 - 180)

    lat_r = np.radians(lat)
    lon_r = np.radians(lon)

    k0 = 0.9996
    ep = ecc**2 / (1 - ecc**2)

    N = a / np.sqrt(1 - ecc**2 * np.sin(lat_r)**2)
    T = np.tan(lat_r)**2
    C = ep * np.cos(lat_r)**2
    A = (lon_r - lon0) * np.cos(lat_r)
    M = a * (
        (1 - ecc**2/4 - 3*ecc**4/64 - 5*ecc**6/256) * lat_r
        - (3*ecc**2/8 + 3*ecc**4/32 + 45*ecc**6/1024) * np.sin(2*lat_r)
        + (15*ecc**4/256 + 45*ecc**6/1024) * np.sin(4*lat_r)
        - 35*ecc**6/3072 * np.sin(6*lat_r)
    )

    easting = k0 * N * (
        A + (1 - T + C) * A**3 / 6
        + (5 - 18*T + T**2 + 72*C - 58*ep) * A**5 / 120
    ) + feast

    northing = k0 * (
        M + N * np.tan(lat_r) * (
            A**2 / 2
            + (5 - T + 9*C + 4*C**2) * A**4 / 24
            + (61 - 58*T + T**2 + 600*C - 300*ep) * A**6 / 720
        )
    ) + fnorth

    zone = np.where(lat < 0, -zone, zone)
    return easting, northing, zone


# ---------------------------------------------------------------------------
# utm2ll
# ---------------------------------------------------------------------------

def utm2ll(
    easting:  float | np.ndarray,
    northing: float | np.ndarray,
    zone:     int   | np.ndarray,
    datum:    str = 'wgs84',
) -> tuple[np.ndarray, np.ndarray]:
    """Convert UTM coordinates to geographic (lon, lat).

    Parameters
    ----------
    easting, northing : array-like
        UTM coordinates in metres.
    zone : int or array-like
        UTM zone number.  Negative values indicate the southern hemisphere.
    datum : str
        Reference ellipsoid name (default ``'wgs84'``).

    Returns
    -------
    lon, lat : ndarray
        Decimal degrees.
    """
    easting  = np.asarray(easting,  dtype=float)
    northing = np.asarray(northing, dtype=float)
    zone     = np.asarray(zone,     dtype=int)

    zabs = np.abs(zone)

    a, b = _get_datum(datum)
    ecc  = _eccentricity(a, b)
    feast, fnorth = _false_origin(datum)

    lon0 = np.radians(6 * (zabs - 1) + 3 - 180)

    k0 = 0.9996
    ep = ecc**2 / (1 - ecc**2)
    e1 = (1 - np.sqrt(1 - ecc**2)) / (1 + np.sqrt(1 - ecc**2))

    x = easting  - feast
    y = northing - fnorth

    M   = y / k0
    mu  = M / (a * (1 - ecc**2/4 - 3*ecc**4/64 - 5*ecc**6/256 - 7*ecc**8/1024))

    lat1 = (
        mu
        + (3*e1/2   - 27*e1**3/32)  * np.sin(2*mu)
        + (21*e1**2/16 - 55*e1**4/32) * np.sin(4*mu)
        + (151*e1**3/96)              * np.sin(6*mu)
        + (1097*e1**4/512)            * np.sin(8*mu)
    )

    C1 = ep * np.cos(lat1)**2
    T1 = np.tan(lat1)**2
    N1 = a / np.sqrt(1 - ecc**2 * np.sin(lat1)**2)
    R1 = a * (1 - ecc**2) / (1 - ecc**2 * np.sin(lat1)**2)**1.5
    D  = x / (N1 * k0)

    lat = np.degrees(
        lat1
        - (N1 * np.tan(lat1) / R1) * (
            D**2 / 2
            - (5 + 3*T1 + 10*C1 - 4*C1**2 - 9*ep) * D**4 / 24
            + (61 + 90*T1 + 298*C1 + 45*T1**2 - 252*ep - 3*C1**2) * D**6 / 720
        )
    )

    lon = np.degrees(
        lon0
        + (D - (1 + 2*T1 + C1) * D**3 / 6
           + (5 - 2*C1 + 28*T1 - 3*C1**2 + 8*ep + 24*T1**2) * D**5 / 120)
        / np.cos(lat1)
    )

    return lon, lat


# ---------------------------------------------------------------------------
# Spherical geometry utilities
# ---------------------------------------------------------------------------

def sphangle(P: tuple, Q: tuple) -> float:
    """Great-circle angular distance between P and Q (degrees).

    Parameters
    ----------
    P, Q : (lon, lat) in decimal degrees.
    """
    lon1, lat1 = np.radians(P[0]), np.radians(P[1])
    lon2, lat2 = np.radians(Q[0]), np.radians(Q[1])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat/2)**2 + np.cos(lat1)*np.cos(lat2)*np.sin(dlon/2)**2
    return np.degrees(2 * np.arcsin(np.sqrt(a)))


def _rotmat(lon_deg: float, lat_deg: float, angle_deg: float) -> np.ndarray:
    """3-D rotation matrix: rotate *angle_deg* about the axis at (lon, lat)."""
    lon = np.radians(lon_deg)
    lat = np.radians(lat_deg)
    a   = np.radians(angle_deg)

    # Euler axis in Cartesian
    u = np.array([
        np.cos(lat) * np.cos(lon),
        np.cos(lat) * np.sin(lon),
        np.sin(lat),
    ])
    ca, sa = np.cos(a), np.sin(a)
    ux, uy, uz = u

    return np.array([
        [ca + ux**2*(1-ca),    ux*uy*(1-ca) - uz*sa, ux*uz*(1-ca) + uy*sa],
        [uy*ux*(1-ca) + uz*sa, ca + uy**2*(1-ca),    uy*uz*(1-ca) - ux*sa],
        [uz*ux*(1-ca) - uy*sa, uz*uy*(1-ca) + ux*sa, ca + uz**2*(1-ca)  ],
    ])


def gcpoints(
    P: tuple,
    Q: tuple,
    n: int = 100,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute *n* points along the great-circle arc from P to Q.

    Parameters
    ----------
    P, Q : (lon, lat) in decimal degrees.
    n : int
        Number of output points (including endpoints).

    Returns
    -------
    lon, lat : ndarray
        Coordinates of the *n* arc points in decimal degrees.
    R : ndarray, shape (3, 3)
        Combined rotation matrix that maps P → (0, 0) and Q → equator.
    """
    def to_cart(lon_deg, lat_deg):
        lo, la = np.radians(lon_deg), np.radians(lat_deg)
        return np.array([np.cos(la)*np.cos(lo),
                         np.cos(la)*np.sin(lo),
                         np.sin(la)])

    def to_sph(v):
        r = np.linalg.norm(v, axis=0)
        lat = np.degrees(np.arcsin(v[2] / r))
        lon = np.degrees(np.arctan2(v[1], v[0]))
        return lon, lat

    # Step 1: Rotate P to (0, 0)
    R1 = _rotmat(0, -90, P[0])   # rotate about south pole by lon(P)
    R2 = _rotmat(90,  0, P[1])   # rotate about (90°, 0) by lat(P)

    q  = to_cart(*Q)
    qp = R2 @ R1 @ q
    Qp_lon, Qp_lat = to_sph(qp)

    # Step 2: Rotate Q' to the equator
    with np.errstate(divide='ignore', invalid='ignore'):
        cosQp = np.cos(np.radians(Qp_lon)) * np.cos(np.radians(Qp_lat))
        ang   = np.arccos(cosQp)
        Omega = np.degrees(
            np.arccos(np.tan(np.radians(Qp_lon)) / np.tan(ang))
        )
    if Qp_lat > 0:
        Omega = -Omega

    R3 = _rotmat(0, 0, Omega)
    R  = R3 @ R2 @ R1

    # Q'' longitude on equator
    qpp = R @ q
    Qpp_lon, _ = to_sph(qpp)

    # Step 3: Generate n equatorial points (0,0) → (Qpp_lon, 0) and un-rotate
    lons_eq = np.linspace(0, np.radians(Qpp_lon), n)
    eq_pts  = np.row_stack([
        np.cos(lons_eq),
        np.sin(lons_eq),
        np.zeros(n),
    ])

    arc = R.T @ eq_pts  # un-rotate
    lon, lat = to_sph(arc)
    return lon, lat, R


# ---------------------------------------------------------------------------
# Authalic latitude
# ---------------------------------------------------------------------------

def _q_authalic(ecc: float, phi: np.ndarray) -> np.ndarray:
    """Snyder q function (eq. 3-12) used in authalic latitude computation."""
    s = np.sin(phi)
    with np.errstate(divide='ignore', invalid='ignore'):
        return (1 - ecc**2) * (
            s / (1 - ecc**2 * s**2)
            - np.log((1 - ecc * s) / (1 + ecc * s)) / (2 * ecc)
        )


def authalic_lat(
    lat: float | np.ndarray,
    datum: str = 'wgs84',
) -> tuple[np.ndarray, float]:
    """Convert geodetic latitude to authalic latitude.

    The authalic sphere is the equal-area sphere approximating the reference
    ellipsoid (Snyder 1987, pp. 16–17, eqs 3-11 to 3-13).

    Parameters
    ----------
    lat : array-like
        Geodetic latitude in decimal degrees.
    datum : str
        Reference ellipsoid (default ``'wgs84'``).

    Returns
    -------
    alat : ndarray
        Authalic latitude in decimal degrees.
    Re : float
        Radius of the equal-area sphere in metres.
    """
    lat = np.asarray(lat, dtype=float)
    a, b = _get_datum(datum)
    ecc  = _eccentricity(a, b)
    phi  = np.radians(lat)
    q    = _q_authalic(ecc, phi)
    qp   = float(_q_authalic(ecc, np.pi / 2))
    alat = np.degrees(np.arcsin(np.clip(q / qp, -1.0, 1.0)))
    Re   = a * np.sqrt(qp / 2)
    return alat, Re


def _geodetic_from_authalic(beta: np.ndarray, ecc: float) -> np.ndarray:
    """Series: authalic latitude (radians) → geodetic latitude (radians).

    Accurate to < 0.5 mm for WGS84 (Bowring 1983).
    """
    e2, e4, e6 = ecc**2, ecc**4, ecc**6
    return (
        beta
        + (e2/3      + 31*e4/180   + 517*e6/5040) * np.sin(2*beta)
        + (23*e4/360 + 251*e6/3780)               * np.sin(4*beta)
        + (761*e6/45360)                           * np.sin(6*beta)
    )


# ---------------------------------------------------------------------------
# Mercator projection
# ---------------------------------------------------------------------------

def mercator(
    lon: float | np.ndarray,
    lat: float | np.ndarray,
    datum: str = 'wgs84',
    lon0: float = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Forward ellipsoidal Mercator projection (authalic-sphere form).

    Parameters
    ----------
    lon, lat : array-like
        Geographic coordinates in decimal degrees.
    datum : str
        Reference ellipsoid (default ``'wgs84'``).
    lon0 : float
        Central meridian in decimal degrees (default 0).

    Returns
    -------
    x, y : ndarray
        Map coordinates in metres.
    """
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    alat, Re = authalic_lat(lat, datum)
    alat_r = np.radians(alat)
    x = Re * np.radians(lon - lon0)
    with np.errstate(divide='ignore', invalid='ignore'):
        y = Re * np.log(np.tan(np.pi / 4 + alat_r / 2))
    return x, y


def inv_mercator(
    x: float | np.ndarray,
    y: float | np.ndarray,
    datum: str = 'wgs84',
    lon0: float = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Inverse ellipsoidal Mercator projection (authalic-sphere form).

    Parameters
    ----------
    x, y : array-like
        Map coordinates in metres.
    datum : str
        Reference ellipsoid (default ``'wgs84'``).
    lon0 : float
        Central meridian in decimal degrees (default 0).

    Returns
    -------
    lon, lat : ndarray
        Geographic coordinates in decimal degrees.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    a, b = _get_datum(datum)
    ecc  = _eccentricity(a, b)
    Re   = a * np.sqrt(float(_q_authalic(ecc, np.pi / 2)) / 2)
    lon  = np.degrees(x / Re) + lon0
    alat = 2 * np.arctan(np.exp(y / Re)) - np.pi / 2
    lat  = np.degrees(_geodetic_from_authalic(alat, ecc))
    return lon, lat


# ---------------------------------------------------------------------------
# Polar stereographic projection (NSIDC / SCAR convention)
# ---------------------------------------------------------------------------

_WGS84_A = 6_378_137.0    # semi-major axis [m]
_WGS84_E = 0.08181919     # eccentricity


def polarstereo_fwd(
    lon:   float | np.ndarray,
    lat:   float | np.ndarray,
    lat_c: float = -70.0,
    lon0:  float =   0.0,
    a:     float = _WGS84_A,
    e:     float = _WGS84_E,
) -> tuple[np.ndarray, np.ndarray]:
    """Polar stereographic forward projection.

    Used by NSIDC and SCAR for Antarctic/Arctic sea-ice grids.

    Parameters
    ----------
    lon, lat : array-like
        Geographic coordinates in decimal degrees.
    lat_c : float
        Standard parallel (latitude of true scale) in degrees.
        Negative for south-polar aspect (default ``-70``).
    lon0 : float
        Central meridian (positive-Y axis) in degrees (default ``0``).
    a : float
        Semi-major axis in metres (default WGS84).
    e : float
        Eccentricity (default WGS84).

    Returns
    -------
    x, y : ndarray
        Map coordinates in metres.

    References
    ----------
    Snyder (1987) pp. 154–163, polar stereographic with known phi_c.
    """
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)

    phi   = np.radians(lat)
    lam   = np.radians(lon)
    phi_c = np.radians(lat_c)
    lam0  = np.radians(lon0)

    if phi_c < 0:
        pm = -1
        phi, phi_c, lam, lam0 = -phi, -phi_c, -lam, -lam0
    else:
        pm = 1

    def _t(p):
        return np.tan(np.pi / 4 - p / 2) / (
            (1 - e * np.sin(p)) / (1 + e * np.sin(p))
        ) ** (e / 2)

    t_c = float(_t(phi_c))
    m_c = float(np.cos(phi_c) / np.sqrt(1 - e**2 * np.sin(phi_c)**2))
    rho = a * m_c * _t(phi) / t_c

    x = pm * rho * np.sin(lam - lam0)
    y = -pm * rho * np.cos(lam - lam0)
    return x, y


def polarstereo_inv(
    x:     float | np.ndarray,
    y:     float | np.ndarray,
    lat_c: float = -70.0,
    lon0:  float =   0.0,
    a:     float = _WGS84_A,
    e:     float = _WGS84_E,
) -> tuple[np.ndarray, np.ndarray]:
    """Polar stereographic inverse projection.

    Parameters
    ----------
    x, y : array-like
        Map coordinates in metres.
    lat_c : float
        Standard parallel in degrees (default ``-70``).
    lon0 : float
        Central meridian in degrees (default ``0``).
    a : float
        Semi-major axis in metres (default WGS84).
    e : float
        Eccentricity (default WGS84).

    Returns
    -------
    lon, lat : ndarray
        Geographic coordinates in decimal degrees.

    References
    ----------
    Snyder (1987) pp. 154–163.  Series form (not iterative).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    phi_c = np.radians(lat_c)
    lam0  = np.radians(lon0)

    if phi_c < 0:
        pm = -1
        phi_c, lam0, x, y = -phi_c, -lam0, -x.copy(), -y.copy()
    else:
        pm = 1

    t_c = np.tan(np.pi / 4 - phi_c / 2) / (
        (1 - e * np.sin(phi_c)) / (1 + e * np.sin(phi_c))
    ) ** (e / 2)
    m_c  = np.cos(phi_c) / np.sqrt(1 - e**2 * np.sin(phi_c)**2)
    rho  = np.sqrt(x**2 + y**2)
    t    = rho * t_c / (a * m_c)

    # Series expansion (Snyder p. 162) instead of iteration
    chi  = np.pi / 2 - 2 * np.arctan(t)
    e2, e4, e6, e8 = e**2, e**4, e**6, e**8
    phi = (
        chi
        + (e2/2   + 5*e4/24   + e6/12     + 13*e8/360)  * np.sin(2*chi)
        + (7*e4/48  + 29*e6/240 + 811*e8/11520)          * np.sin(4*chi)
        + (7*e6/120 + 81*e8/1120)                        * np.sin(6*chi)
        + (4279*e8/161280)                               * np.sin(8*chi)
    )
    lam  = lam0 + np.arctan2(x, -y)

    phi = pm * phi
    lam = pm * lam
    lam = np.mod(lam + np.pi, 2 * np.pi) - np.pi

    return np.degrees(lam), np.degrees(phi)


# ---------------------------------------------------------------------------
# Small circles
# ---------------------------------------------------------------------------

def smallcircle(
    lon0:  float,
    lat0:  float,
    delta: float,
    n:     int = 360,
) -> tuple[np.ndarray, np.ndarray]:
    """Points on the small circle of angular radius *delta* about (lon0, lat0).

    Parameters
    ----------
    lon0, lat0 : float
        Centre of the small circle in decimal degrees.
    delta : float
        Angular radius in degrees.
    n : int
        Number of output points (default 360).

    Returns
    -------
    clon, clat : ndarray
        Longitudes and latitudes in decimal degrees.
    """
    # Build circle around north pole, then rotate to (lon0, lat0)
    clon_r = np.linspace(-np.pi, np.pi, n, endpoint=False)
    clat_r = np.radians(90.0 - delta) * np.ones(n)

    q = np.row_stack([
        np.cos(clat_r) * np.cos(clon_r),
        np.cos(clat_r) * np.sin(clon_r),
        np.sin(clat_r),
    ])

    # R maps (lon0, lat0) → north pole; R.T is the inverse
    R  = _rotmat(lon0 - 90, 0, 90 - lat0)
    qp = R.T @ q

    clon = np.degrees(np.arctan2(qp[1], qp[0]))
    clat = np.degrees(np.arcsin(np.clip(qp[2], -1.0, 1.0)))
    return clon, clat


# ---------------------------------------------------------------------------
# Area calculations
# ---------------------------------------------------------------------------

def pixarea(
    lat1: float | np.ndarray,
    lat2: float | np.ndarray,
    lon1: float | np.ndarray,
    lon2: float | np.ndarray,
    datum: str | None = None,
) -> np.ndarray:
    """Fractional area of a lat/lon quad pixel on a sphere.

    Returns the area as a dimensionless fraction of the total sphere surface
    (integrates to 1 over the entire globe).  Multiply by 4π R² to convert
    to absolute area.

    Parameters
    ----------
    lat1, lat2 : array-like
        Southern and northern bounding latitudes in degrees.
    lon1, lon2 : array-like
        Western and eastern bounding longitudes in degrees.
    datum : str or None
        When given, latitudes are converted to authalic latitude before
        the area calculation (ellipsoidal correction).  Use ``None``
        (default) for the spherical formula.

    Returns
    -------
    A : ndarray
        Dimensionless fractional area.
    """
    lat1 = np.asarray(lat1, dtype=float)
    lat2 = np.asarray(lat2, dtype=float)
    lon1 = np.asarray(lon1, dtype=float)
    lon2 = np.asarray(lon2, dtype=float)

    if datum is not None:
        lat1, _ = authalic_lat(lat1, datum)
        lat2, _ = authalic_lat(lat2, datum)

    r1   = np.radians(lat1)
    r2   = np.radians(lat2)
    dlon = np.radians(np.abs(lon2 - lon1))
    return dlon / (4 * np.pi) * np.abs(np.sin(r2) - np.sin(r1))


def rectarea(
    lonlim: tuple | np.ndarray,
    latlim: tuple | np.ndarray,
    R: float = 6371.0,
) -> float:
    """Area of a lat/lon rectangle on a sphere.

    Parameters
    ----------
    lonlim : (lon_min, lon_max)
        Longitude limits in degrees.
    latlim : (lat_min, lat_max)
        Latitude limits in degrees.
    R : float
        Sphere radius in km (default 6371 km).

    Returns
    -------
    float
        Area in km².
    """
    lat1 = np.radians(min(latlim))
    lat2 = np.radians(max(latlim))
    dlon = np.radians(abs(lonlim[1] - lonlim[0]))
    return R**2 * dlon * abs(np.sin(lat2) - np.sin(lat1))
