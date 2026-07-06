"""
ECM1 global crustal structure model.

Provides 1°×1° global estimates of crustal layer thicknesses, seismic
velocities, and densities.

Typical usage
-------------
>>> from src.geophys.ecm1 import load_ecm1, interpolate_to_points, plot_ecm1
>>> df = load_ecm1()
>>> hc = interpolate_to_points(lons, lats, 'Hc')   # crustal thickness (km)
>>> ax = plot_ecm1('Hc', projection='robinson')
"""

from __future__ import annotations

import functools
import json
import pathlib

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Column metadata
# ---------------------------------------------------------------------------

COLUMN_INFO: dict[str, str] = {
    'Hcc':  'Crystalline crust thickness (km)',
    'Sed':  'Sediment thickness (km)',
    'Hc':   'Total crust thickness (km)',
    'Type': 'Crust type code',
    'DLy1': 'Upper crust layer thickness (km)',
    'DLy2': 'Middle crust layer thickness (km)',
    'DLy3': 'Lower crust layer thickness (km)',
    'TLy1': 'Upper crust base temperature (°C)',
    'TLy2': 'Middle crust base temperature (°C)',
    'TLy3': 'Lower crust base temperature (°C)',
    'VP1':  'Upper crust P-wave velocity (km/s)',
    'VP2':  'Middle crust P-wave velocity (km/s)',
    'VP3':  'Lower crust P-wave velocity (km/s)',
    'VS1':  'Upper crust S-wave velocity (km/s)',
    'VS2':  'Middle crust S-wave velocity (km/s)',
    'VS3':  'Lower crust S-wave velocity (km/s)',
    'VPN':  'Pn velocity — uppermost mantle (km/s)',
    'VSN':  'Sn velocity — uppermost mantle (km/s)',
    'RHO1': 'Upper crust density (g/cm³)',
    'RHO2': 'Middle crust density (g/cm³)',
    'RHO3': 'Lower crust density (g/cm³)',
    'RHON': 'Upper mantle density (g/cm³)',
}

_COLUMNS = [
    'Numb', 'Lon', 'Lat', 'Hcc', 'Sed', 'Hc', 'Type',
    'DLy1', 'DLy2', 'DLy3', 'TLy1', 'TLy2', 'TLy3',
    'VP1', 'VP2', 'VP3', 'VS1', 'VS2', 'VS3', 'VPN', 'VSN',
    'RHO1', 'RHO2', 'RHO3', 'RHON',
]

# Regular 1°×1° grid cell centres
_LONS     = np.arange(-179.5, 180.0, 1.0)    # 360 values, west→east
_LATS_ASC = np.arange(-89.5,   90.0, 1.0)    # 180 values, south→north
#   File stores data north→south; we flip on load so arrays index [lat_asc, lon]


# ---------------------------------------------------------------------------
# Path resolution
# ---------------------------------------------------------------------------

def _resolve_ecm1_path() -> pathlib.Path:
    """Find ECM1.txt via data_registry.json or raise FileNotFoundError."""
    registry_path = (
        pathlib.Path(__file__).resolve().parent.parent.parent
        / 'data_access' / 'data_registry.json'
    )
    if registry_path.exists():
        with registry_path.open() as f:
            registry = json.load(f)
        for ds in registry.get('datasets', []):
            if ds.get('name') == 'ECM1':
                loc = pathlib.Path(ds['location']).expanduser()
                if loc.exists():
                    return loc
    raise FileNotFoundError(
        'ECM1.txt not found. Update the "location" field for "ECM1" in '
        'data_access/data_registry.json or pass path= explicitly to load_ecm1().'
    )


# ---------------------------------------------------------------------------
# Core loader (cached)
# ---------------------------------------------------------------------------

@functools.lru_cache(maxsize=1)
def load_ecm1(path: str | None = None) -> pd.DataFrame:
    """Load the ECM1 crustal model into a DataFrame.

    The result is cached after the first call; subsequent calls with the
    same *path* are instant.

    Parameters
    ----------
    path : str, optional
        Path to ECM1.txt.  Resolved from data_registry.json when *None*.

    Returns
    -------
    pandas.DataFrame
        64,800 rows × 25 columns.  See :data:`COLUMN_INFO` for descriptions.
    """
    p = pathlib.Path(path) if path is not None else _resolve_ecm1_path()
    df = pd.read_csv(p, sep=r'\s+', names=_COLUMNS, header=0)
    return df


# ---------------------------------------------------------------------------
# Grid extraction
# ---------------------------------------------------------------------------

def to_grid(column: str, path: str | None = None):
    """Return an ECM1 field as a 2-D grid for :func:`~matplotlib.axes.Axes.pcolormesh`.

    Parameters
    ----------
    column : str
        ECM1 column name (e.g. ``'Hc'`` for total crustal thickness).
    path : str, optional
        Passed to :func:`load_ecm1`.

    Returns
    -------
    lons : (361,) ndarray
        Cell-edge longitudes (−180 → 180).
    lats : (181,) ndarray
        Cell-edge latitudes (−90 → 90).
    grid : (180, 360) ndarray
        Field values on the 1°×1° grid, indexed ``[lat, lon]`` with lat
        increasing south→north.
    """
    df = load_ecm1(path)
    grid = df[column].values.reshape(180, 360)   # file order: lat descending
    grid = grid[::-1, :]                          # flip to lat ascending
    lons = np.arange(-180.0, 181.0, 1.0)          # 361 cell edges
    lats = np.arange(-90.0,   91.0, 1.0)          # 181 cell edges
    return lons, lats, grid


# ---------------------------------------------------------------------------
# Interpolation
# ---------------------------------------------------------------------------

def interpolate_to_points(
    lon,
    lat,
    column: str,
    method: str = 'linear',
    path: str | None = None,
) -> np.ndarray:
    """Interpolate an ECM1 field to arbitrary geographic points.

    Parameters
    ----------
    lon, lat : array-like
        Query longitudes / latitudes in decimal degrees.  NaN values are
        propagated as NaN in the output.
    column : str
        ECM1 column to interpolate (e.g. ``'Hc'``).
    method : {'linear', 'nearest'}
        Interpolation method passed to
        :class:`~scipy.interpolate.RegularGridInterpolator`.
    path : str, optional
        Passed to :func:`load_ecm1`.

    Returns
    -------
    numpy.ndarray
        Interpolated values, same length as the input arrays.
    """
    from scipy.interpolate import RegularGridInterpolator

    df = load_ecm1(path)
    grid = df[column].values.reshape(180, 360)[::-1, :]  # lat ascending

    interp = RegularGridInterpolator(
        (_LATS_ASC, _LONS), grid,
        method=method,
        bounds_error=False,
        fill_value=np.nan,
    )

    lon = np.asarray(lon, dtype=float).ravel()
    lat = np.asarray(lat, dtype=float).ravel()

    # Clamp longitudes to [-180, 180] to handle slight out-of-range values
    lon = ((lon + 180) % 360) - 180

    valid = np.isfinite(lon) & np.isfinite(lat)
    result = np.full(len(lon), np.nan)
    if valid.any():
        result[valid] = interp(np.column_stack([lat[valid], lon[valid]]))
    return result


# ---------------------------------------------------------------------------
# Mapping
# ---------------------------------------------------------------------------

def plot_ecm1(
    column: str = 'Hc',
    projection: str = 'geographic',
    lon0: float = 0.0,
    lat0: float = 0.0,
    ax=None,
    colorbar: bool = True,
    cmap: str = 'viridis',
    path: str | None = None,
    **kw,
):
    """Plot an ECM1 field as a pseudocolor map.

    Overlays correctly on axes created by :func:`~src.maptools.coast.plotcoast`
    using the same *projection*, *lon0*, and *lat0*.

    Parameters
    ----------
    column : str
        ECM1 column to plot.  Default ``'Hc'`` (total crustal thickness).
    projection : str
        Projection name — same options as :func:`~src.maptools.coast.plotcoast`.
    lon0, lat0 : float
        Map centre (must match the base map).
    ax : matplotlib.axes.Axes, optional
        Target axes.  Created if *None*.
    colorbar : bool
        Attach a colorbar to *ax*.
    cmap : str
        Matplotlib colormap name.
    path : str, optional
        Passed to :func:`load_ecm1`.
    **kw
        Forwarded to :func:`~matplotlib.axes.Axes.pcolormesh`
        (e.g. ``vmin``, ``vmax``, ``shading``).

    Returns
    -------
    matplotlib.axes.Axes
    """
    import matplotlib.pyplot as plt

    # Lazy import so the module loads without requiring maptools
    import sys as _sys
    import pathlib as _pathlib
    _root = _pathlib.Path(__file__).resolve().parent.parent.parent
    if str(_root) not in _sys.path:
        _sys.path.insert(0, str(_root))
    from src.maptools.coast import _project

    if ax is None:
        _, ax = plt.subplots()

    kw.setdefault('shading', 'flat')

    lons_e, lats_e, grid = to_grid(column, path)

    if projection == 'geographic':
        # Shift so lon0 is at centre
        lon_mesh, lat_mesh = np.meshgrid(lons_e - lon0, lats_e)
        pcm = ax.pcolormesh(lon_mesh, lat_mesh, grid, cmap=cmap, **kw)
    else:
        # Project the (N+1)×(M+1) edge mesh; pcolormesh handles irregular quads
        lon_mesh, lat_mesh = np.meshgrid(lons_e, lats_e)
        x, y = _project(
            lon_mesh.ravel().astype(float),
            lat_mesh.ravel().astype(float),
            projection, lon0, lat0,
        )
        x = x.reshape(lon_mesh.shape)
        y = y.reshape(lat_mesh.shape)
        # Mask cells that project to NaN (outside map boundary)
        valid_mask = np.isfinite(x[:-1, :-1]) & np.isfinite(y[:-1, :-1])
        masked_grid = np.where(valid_mask, grid, np.nan)
        pcm = ax.pcolormesh(x, y, masked_grid, cmap=cmap, **kw)

    if colorbar:
        label = COLUMN_INFO.get(column, column)
        plt.colorbar(pcm, ax=ax, label=label, shrink=0.7, pad=0.02)

    return ax
