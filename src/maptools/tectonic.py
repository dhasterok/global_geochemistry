"""
Tectonic and geographic boundary data from global_tectonics.

Loads shapefiles from the local clone of
https://github.com/dhasterok/global_tectonics
and projects/plots them on map axes produced by :func:`~coast.plotcoast`.

Registered datasets
-------------------
``'plate_boundaries'``
    Tectonic plate boundary polylines (spreading centres, subduction zones,
    transform faults, etc.).  571 features.  Field ``'type'`` values:
    ``'subduction zone'``, ``'spreading center'``, ``'collision zone'``,
    ``'dextral transform'``, ``'sinistral transform'``, ``'extension zone'``,
    ``'inferred'``.

``'plates'``
    Tectonic plate polygons.  280 features.

``'oc_boundaries'``
    Ocean–continent boundary polylines.  93 features.

``'provinces'``
    Global geologic provinces (polygons).  914 features.

``'cratons'``
    Craton polygons.  114 features.

``'coastline_hi'``
    High-resolution GSHHS coastline polygons.  6042 features.

Usage
-----
>>> from global_geochemistry.maptools.coast import plotcoast
>>> from global_geochemistry.maptools.tectonic import plotboundaries, plotpolygons
>>> ax = plotcoast(projection='robinson')
>>> plotpolygons('plates', projection='robinson', ax=ax, zorder=0)
>>> plotboundaries('plate_boundaries', projection='robinson', ax=ax,
...                attribute='type', value='subduction zone',
...                color='k', barbs=True)
"""

from __future__ import annotations

import colorsys
import functools
import json
import pathlib
import warnings

import numpy as np
import matplotlib.pyplot as plt

from shapely.geometry import box as _shp_box, Polygon as _ShpPolygon, MultiPolygon as _ShpMultiPolygon

from .coast import (_preprocess_lines, _project, _clip_lines,
                    _clip_ring_at_folds, _presplit_at_equator,
                    _WATERMAN_INNER_FOLDS,
                    _clip_to_panel, _waterman_panels, _densify_panel_boundary)

# ---------------------------------------------------------------------------
# QGIS plate colours (from global_tectonics/styles/plates_types.qml)
# ---------------------------------------------------------------------------

_QGIS_PLATE_COLORS: dict[str, tuple[float, float, float]] = {
    name: tuple(c / 255 for c in rgb)
    for name, rgb in {
        'African Plate':        (229, 216, 202),
        'Antarctic Plate':      (245, 218, 233),
        'South American Plate': (251, 209, 206),
        'Arabian Plate':        (204, 224, 237),
        'Australian Plate':     (255, 229, 206),
        'Caribbean Plate':      (255, 239, 207),
        'Cocos Plate':          (228, 220, 237),
        'Eurasian Plate':       (254, 248, 223),
        'Indian Plate':         (228, 201, 200),
        'Juan de Fuca Plate':   (245, 224, 205),
        'North American Plate': (208, 232, 206),
        'Nazca Plate':          (207, 236, 239),
        'Pacific Plate':        (209, 220, 241),
        'Philippine Plate':     (246, 248, 222),
        'Scotia Plate':         (210, 219, 224),
        'Somali Plate':         (237, 237, 237),
    }.items()
}


# ---------------------------------------------------------------------------
# Data registry — path resolved from data_registry.json at runtime
# ---------------------------------------------------------------------------

def _get_tectonics_root() -> pathlib.Path:
    """Return the global_tectonics data root, resolved from the data registry.

    Search order:
    1. ``"location"`` field in data_registry.json for the ``"global_tectonics"``
       dataset (set automatically by :func:`~data_access.download_data.download_dataset`).
    2. ``<repo_root>/data/global_tectonics/``
    3. ``~/global_tectonics/``

    Raises
    ------
    FileNotFoundError
        If no valid location is found.  Run
        ``data_access.download_data.download_dataset('global_tectonics')``
        to download the data automatically.
    """
    _repo = pathlib.Path(__file__).resolve().parent.parent.parent
    registry_path = _repo / 'data_access' / 'data_registry.json'

    if registry_path.exists():
        with registry_path.open() as f:
            registry = json.load(f)
        for ds in registry.get('datasets', []):
            if ds.get('name') == 'global_tectonics':
                loc = ds.get('location', '')
                if loc and loc != 'auto':
                    p = pathlib.Path(loc).expanduser()
                    if p.exists():
                        return p

    # Well-known fallback locations
    for candidate in [
        _repo / 'data' / 'global_tectonics',
        pathlib.Path.home() / 'global_tectonics',
    ]:
        if candidate.exists():
            return candidate

    raise FileNotFoundError(
        'global_tectonics data not found.\n'
        'Download it with:\n'
        '    from data_access.download_data import download_dataset\n'
        '    download_dataset("global_tectonics")\n'
        'or set its "location" in data_access/data_registry.json.'
    )


def _build_registry() -> dict[str, pathlib.Path]:
    root = _get_tectonics_root()
    return {
        'plate_boundaries': root / 'plates&provinces/shp/plate_boundaries.shp',
        'plates':           root / 'plates&provinces/shp/plates.shp',
        'oc_boundaries':    root / 'plates&provinces/shp/oc_boundaries.shp',
        'provinces':        root / 'plates&provinces/shp/global_gprv.shp',
        'cratons':          root / 'plates&provinces/shp/cratons.shp',
        'coastline_hi':     root / 'polygon_data/GSHHS_I_L1.shp',
    }


# Module-level cache — populated on first use
_REGISTRY: dict[str, pathlib.Path] = {}


def _get_registry() -> dict[str, pathlib.Path]:
    global _REGISTRY
    if not _REGISTRY:
        _REGISTRY = _build_registry()
    return _REGISTRY


def list_datasets() -> list[str]:
    """Return names of all registered tectonic datasets."""
    return sorted(_get_registry())


# ---------------------------------------------------------------------------
# Shapefile I/O helpers
# ---------------------------------------------------------------------------

def _open_shapefile(name_or_path):
    try:
        import shapefile as pyshp
    except ImportError:
        raise ImportError(
            'pyshp is required for shapefile loading.  '
            'Install it with:  pip install pyshp'
        ) from None
    path = _get_registry().get(str(name_or_path), pathlib.Path(name_or_path))
    if not pathlib.Path(path).exists():
        raise FileNotFoundError(
            f"Shapefile not found: {path}\n"
            f"Registered names: {list_datasets()}"
        )
    return pyshp.Reader(str(path))


def _filter_index(sf, attribute):
    """Return (field_names, filter_idx) for attribute filtering."""
    field_names = [f[0] for f in sf.fields[1:]]
    if attribute is None:
        return field_names, None
    if attribute not in field_names:
        raise ValueError(
            f"Attribute '{attribute}' not found.  "
            f"Available fields: {field_names}"
        )
    return field_names, field_names.index(attribute)


def _split_at_nan(x, y):
    """Split NaN-separated coordinate arrays into finite (x, y) segment lists."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    segments = []
    start = None
    for i in range(len(x)):
        if valid[i]:
            if start is None:
                start = i
        else:
            if start is not None and i - start >= 2:
                segments.append((x[start:i], y[start:i]))
            start = None
    if start is not None and len(x) - start >= 2:
        segments.append((x[start:], y[start:]))
    return segments


# ---------------------------------------------------------------------------
# Shapefile loaders
# ---------------------------------------------------------------------------

def load_shapes(
    name_or_path: str | pathlib.Path,
    attribute: str | None = None,
    value=None,
    filters: dict | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Load a shapefile and return NaN-separated (lon, lat) arrays.

    Handles POLYLINE (type 3) and POLYGON (type 5) shapefiles.  Each part
    within each feature is separated by a NaN.

    Parameters
    ----------
    name_or_path : str or Path
        Either a registered dataset name (see :func:`list_datasets`) or a
        path to a ``.shp`` file.
    attribute : str, optional
        If provided together with *value*, only features where the named
        attribute equals *value* are returned.
    value : optional
        Value to match against *attribute*.
    filters : dict, optional
        Additional ``{attribute: value}`` pairs that must ALL match.
        Combined with *attribute*/*value* when both are given.

    Returns
    -------
    lon, lat : numpy.ndarray
        Coordinate arrays in decimal degrees with NaN separators between
        polyline/polygon parts.
    """
    sf = _open_shapefile(name_or_path)
    if sf.shapeType not in (3, 5, 13, 15):
        raise ValueError(
            f"Unsupported shape type {sf.shapeType} ({sf.shapeTypeName}).  "
            "Only POLYLINE and POLYGON are supported."
        )

    # Build combined filter dict from attribute/value + filters
    combined: dict = {}
    if attribute is not None and value is not None:
        combined[attribute] = value
    if filters:
        combined.update(filters)

    field_names = [f[0] for f in sf.fields[1:]]
    filter_pairs = []
    for attr, val in combined.items():
        if attr not in field_names:
            raise ValueError(
                f"Attribute '{attr}' not found.  "
                f"Available fields: {field_names}"
            )
        filter_pairs.append((field_names.index(attr), val))

    lon_segs: list[np.ndarray] = []
    lat_segs: list[np.ndarray] = []

    for rec, sh in zip(sf.iterRecords(), sf.iterShapes()):
        if any(rec[idx] != val for idx, val in filter_pairs):
            continue
        if sh.shapeType == 0:
            continue

        pts = np.asarray(sh.points)
        parts = list(sh.parts) + [len(pts)]

        for k in range(len(parts) - 1):
            seg = pts[parts[k]:parts[k + 1]]
            if len(seg) < 2:
                continue
            lon_segs.append(seg[:, 0])
            lon_segs.append(np.array([np.nan]))
            lat_segs.append(seg[:, 1])
            lat_segs.append(np.array([np.nan]))

    if not lon_segs:
        return np.array([]), np.array([])

    return np.concatenate(lon_segs), np.concatenate(lat_segs)


def load_polygon_parts(
    name_or_path: str | pathlib.Path,
    attribute: str | None = None,
    value=None,
    group_by: str | None = None,
) -> list[tuple[np.ndarray, np.ndarray]] | dict[str, list]:
    """Load polygon shapefile parts as individual (lon, lat) ring arrays.

    Unlike :func:`load_shapes`, each ring is returned as a separate array
    rather than concatenated with NaN separators.  This is suitable for
    :func:`matplotlib.axes.Axes.fill`.

    Parameters
    ----------
    name_or_path : str or Path
        Registered dataset name or path to a ``.shp`` file.
    attribute, value : optional
        Filter to features where ``attribute == value``.
    group_by : str, optional
        Field name to group features by.  If provided, returns a ``dict``
        mapping group value → list of ``(lon, lat)`` arrays.

    Returns
    -------
    list of (lon, lat) arrays  |  dict mapping group → list of (lon, lat)
    """
    sf = _open_shapefile(name_or_path)
    if sf.shapeType not in (5, 15):
        raise ValueError(
            f"Shape type {sf.shapeType} is not a polygon.  "
            "Use load_shapes() for polylines."
        )

    field_names, filter_idx = _filter_index(sf, attribute)
    _, group_idx = _filter_index(sf, group_by) if group_by else (None, None)

    parts_dict: dict = {}
    parts_list: list = []

    for rec, sh in zip(sf.iterRecords(), sf.iterShapes()):
        if filter_idx is not None and rec[filter_idx] != value:
            continue
        if sh.shapeType == 0:
            continue

        pts = np.asarray(sh.points)
        bounds = list(sh.parts) + [len(pts)]

        group_key = str(rec[group_idx]) if group_idx is not None else None

        for k in range(len(bounds) - 1):
            ring = pts[bounds[k]:bounds[k + 1]]
            if len(ring) < 3:
                continue
            lon_r, lat_r = ring[:, 0], ring[:, 1]
            # In ESRI shapefiles exterior rings are CW (negative signed area)
            # and interior rings / holes are CCW (positive signed area).
            # Skip holes so they are not filled as solid polygons.
            area = 0.5 * np.sum(
                lon_r[:-1] * lat_r[1:] - lon_r[1:] * lat_r[:-1]
            )
            if area >= 0:   # CCW = hole → skip
                continue
            entry = (lon_r, lat_r)
            if group_key is not None:
                parts_dict.setdefault(group_key, []).append(entry)
            else:
                parts_list.append(entry)

    return parts_dict if group_key is not None else parts_list


# ---------------------------------------------------------------------------
# Barb drawing
# ---------------------------------------------------------------------------

def _draw_barbs(
    ax,
    x: np.ndarray,
    y: np.ndarray,
    spacing: float | None = None,
    length: float | None = None,
    side: int = 1,
    color='k',
    zorder: int = 4,
):
    """Draw filled triangular barbs perpendicular to a projected polyline.

    Barbs represent the overriding plate on subduction zone traces.  They
    are placed at regular intervals along the line, with each barb's tip
    pointing perpendicular to the local line direction.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    x, y : numpy.ndarray
        Projected coordinate arrays, NaN-separated between segments.
    spacing : float, optional
        Distance between barb bases in map units.  Auto-scaled to
        ~2.5 % of the data bounding box if *None*.
    length : float, optional
        Length of each barb in map units.  Auto-scaled to ~0.8 % of
        the data bounding box if *None*.
    side : int
        ``1`` (default) — barb points to the **right** of the line's
        travel direction.  ``-1`` — points to the left.
        Which physical side corresponds to the overriding plate depends on
        how the trench lines are oriented in the shapefile.
    color : color
        Fill colour for the triangular barbs.
    zorder : int
        Drawing order.  Default 4 (above most map elements).
    """
    segments = _split_at_nan(x, y)
    if not segments:
        return

    # Auto-scale from the bounding box of all projected data
    all_x = np.concatenate([s[0] for s in segments])
    all_y = np.concatenate([s[1] for s in segments])
    data_scale = max(np.ptp(all_x), np.ptp(all_y))

    if spacing is None:
        spacing = data_scale * 0.025
    if length is None:
        length = data_scale * 0.008

    half_base = length * 0.45  # half-width of triangle base along the line

    for xs, ys in segments:
        if len(xs) < 2:
            continue

        dx = np.diff(xs)
        dy = np.diff(ys)
        seg_len = np.sqrt(dx**2 + dy**2)

        # Cumulative arc length
        cum_len = np.concatenate([[0], np.cumsum(seg_len)])
        total_len = cum_len[-1]

        barb_positions = np.arange(spacing / 2, total_len, spacing)

        for pos in barb_positions:
            # Locate the segment index
            idx = min(np.searchsorted(cum_len[1:], pos), len(dx) - 1)
            if seg_len[idx] == 0:
                continue

            t = (pos - cum_len[idx]) / seg_len[idx]
            bx = xs[idx] + t * dx[idx]
            by = ys[idx] + t * dy[idx]

            # Unit vectors: forward and perpendicular (right)
            fx = dx[idx] / seg_len[idx]
            fy = dy[idx] / seg_len[idx]
            px = fy * side    # right perpendicular: rotate CW
            py = -fx * side

            tip_x = bx + length * px
            tip_y = by + length * py
            b1_x = bx - half_base * fx
            b1_y = by - half_base * fy
            b2_x = bx + half_base * fx
            b2_y = by + half_base * fy

            ax.fill([b1_x, b2_x, tip_x], [b1_y, b2_y, tip_y],
                    color=color, zorder=zorder, linewidth=0)


# ---------------------------------------------------------------------------
# Geographic polygon splitting (shapely-based)
# ---------------------------------------------------------------------------

def _iter_polygons(geom):
    """Yield every Polygon component of a shapely geometry."""
    if geom.is_empty:
        return
    t = geom.geom_type
    if t == 'Polygon':
        yield geom
    elif t in ('MultiPolygon', 'GeometryCollection'):
        for g in geom.geoms:
            yield from _iter_polygons(g)


def _geographic_sub_rings(lon_raw, lat_raw, fold_lons=None, split_equator=False):
    """Split a geographic polygon ring at projection fold meridians.

    Uses shapely to clip the ring in geographic space, returning a list of
    closed (lon, lat) sub-rings that each lie entirely within one projection
    cell.  This prevents the antimeridian / fold-meridian artefacts that
    occur when matplotlib tries to fill an open arc.

    Parameters
    ----------
    lon_raw, lat_raw : array_like
        Ring coordinates in decimal degrees.  May be closed (first == last)
        or open.
    fold_lons : list of float
        Meridians at which the projection folds.  Default ``[-180, 180]``
        splits only at the antimeridian.  Pass ``[-180,-90,0,90,180]`` for
        Waterman.
    split_equator : bool
        If *True*, additionally split each sub-ring at lat = 0.  Required
        for projections (e.g. Waterman) that apply a reflection at the
        equator.

    Returns
    -------
    list of (lon, lat) ndarray pairs, each a closed ring.
    """
    try:
        from shapely.geometry import Polygon, box as shp_box
        from shapely.validation import make_valid
    except ImportError:
        raise ImportError(
            'shapely is required for polygon filling.  '
            'Install it with:  pip install shapely'
        ) from None

    if fold_lons is None:
        fold_lons = [-180, 180]

    lon = np.asarray(lon_raw, dtype=float)
    lat = np.asarray(lat_raw, dtype=float)

    # Ensure closed ring
    if not (lon[-1] == lon[0] and lat[-1] == lat[0]):
        lon = np.concatenate([lon, [lon[0]]])
        lat = np.concatenate([lat, [lat[0]]])

    pts_lon = lon[:-1].copy()   # shapely auto-closes; copy so we can snap
    pts_lat = lat[:-1]

    if len(pts_lon) < 3:
        return []

    sorted_folds = sorted(set(fold_lons))

    # Snap vertices within a small tolerance of each fold meridian to the
    # meridian exactly.  Many polygons are drawn deliberately stopping a
    # fraction of a degree short of ±180° (or Waterman fold lons) to avoid
    # QGIS wrapping artefacts; snapping ensures the fill reaches the map frame.
    _SNAP = 0.02   # degrees — catches ±179.999° without disturbing ±179° features
    for fl in sorted_folds:
        pts_lon[np.abs(pts_lon - fl) < _SNAP] = float(fl)

    # Rebuild the closed ring with snapped vertices
    lon = np.concatenate([pts_lon, [pts_lon[0]]])

    # Unwrap longitudes so the ring is continuous (no antimeridian jumps)
    dlon = np.diff(pts_lon)
    dlon_w = ((dlon + 180) % 360) - 180
    lon_cont = np.concatenate([[pts_lon[0]], pts_lon[0] + np.cumsum(dlon_w)])

    lo_cont = lon_cont.min()
    hi_cont = lon_cont.max()

    # Detect crossings in ORIGINAL (post-snap) coordinates.
    # Normal crossings (~±358°) push lo_cont/hi_cont outside the cell range.
    # Degenerate crossings at exactly ±360° (where dlon_w → 0°) are caught by
    # checking all consecutive pairs including the closing segment.
    all_diffs = np.abs(np.diff(pts_lon))
    raw_crossings = bool(np.any(all_diffs > 180 - 1e-9))

    # Check if the unwrapped ring fits entirely within one projection cell
    fits_one_cell = False
    if not raw_crossings:
        for i in range(len(sorted_folds) - 1):
            if lo_cont >= sorted_folds[i] - 1e-9 and hi_cont <= sorted_folds[i + 1] + 1e-9:
                fits_one_cell = True
                break

    if fits_one_cell:
        if not split_equator or lat.min() >= -1e-9 or lat.max() <= 1e-9:
            return [(lon, lat)]
        return _equator_split_ring(lon, lat)

    # Build shapely polygon in unwrapped lon space
    try:
        poly = Polygon(zip(lon_cont, pts_lat))
        if not poly.is_valid:
            poly = make_valid(poly)
        if poly.is_empty or poly.area == 0:
            return []
    except Exception:
        return []

    # Collect all fold boundaries in unwrapped space that intersect [lo_cont, hi_cont]
    boundaries = {lo_cont, hi_cont}
    for fl in sorted_folds:
        k_min = int(np.floor((lo_cont - fl) / 360)) - 1
        k_max = int(np.ceil((hi_cont - fl) / 360)) + 1
        for k in range(k_min, k_max + 1):
            b = fl + k * 360
            if lo_cont - 1e-9 <= b <= hi_cont + 1e-9:
                boundaries.add(b)

    boundaries = sorted(boundaries)
    result = []

    for i in range(len(boundaries) - 1):
        cell_lo = boundaries[i]
        cell_hi = boundaries[i + 1]
        if cell_hi - cell_lo < 1e-9:
            continue

        try:
            clipped = poly.intersection(shp_box(cell_lo - 1e-9, -91,
                                                 cell_hi + 1e-9,  91))
        except Exception:
            continue

        if clipped.is_empty:
            continue

        # Shift clipped coords back to canonical [-180, 180]
        canonical_lo = ((cell_lo + 180) % 360) - 180
        x_shift = canonical_lo - cell_lo

        for sub in _iter_polygons(clipped):
            if sub.is_empty:
                continue
            coords = np.array(sub.exterior.coords)
            sub_lon = coords[:, 0] + x_shift
            sub_lat = coords[:, 1]
            if split_equator:
                result.extend(_equator_split_ring(sub_lon, sub_lat))
            else:
                result.append((sub_lon, sub_lat))

    return result if result else [(lon, lat)]


def _equator_split_ring(lon, lat):
    """Split a ring at lat = 0 for projections with an equatorial reflection."""
    try:
        from shapely.geometry import Polygon, box as shp_box
        from shapely.validation import make_valid
    except ImportError:
        return [(lon, lat)]

    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)

    if lat.min() >= -1e-9 or lat.max() <= 1e-9:
        return [(lon, lat)]

    pts_lon = lon[:-1] if (lon[-1] == lon[0] and lat[-1] == lat[0]) else lon
    pts_lat = lat[:-1] if (lon[-1] == lon[0] and lat[-1] == lat[0]) else lat

    if len(pts_lon) < 3:
        return []

    try:
        poly = Polygon(zip(pts_lon, pts_lat))
        if not poly.is_valid:
            poly = make_valid(poly)
        if poly.is_empty:
            return []
    except Exception:
        return []

    result = []
    for clip_box in [shp_box(-360, 0, 360, 90), shp_box(-360, -90, 360, 0)]:
        try:
            clipped = poly.intersection(clip_box)
        except Exception:
            continue
        for sub in _iter_polygons(clipped):
            if sub.is_empty:
                continue
            coords = np.array(sub.exterior.coords)
            result.append((coords[:, 0], coords[:, 1]))

    return result if result else [(lon, lat)]


# ---------------------------------------------------------------------------
# Tectonic attribute assignment (point-in-polygon)
# ---------------------------------------------------------------------------

# Fields to extract from each dataset
_PLATE_FIELDS    = ['plate', 'plate_code', 'subplate', 'plate_type', 'crust_type']
_PROVINCE_FIELDS = ['prov_name', 'prov_type', 'prov_group', 'lastorogen', 'continent']


@functools.lru_cache(maxsize=4)
def _build_spatial_index(dataset: str):
    """Build and cache a shapely STRtree for a tectonic polygon dataset.

    Returns
    -------
    tree : shapely.STRtree
    geoms : list of shapely geometries (same order as tree)
    records : list of dicts (attribute dicts, same order as tree)
    fields : list of str (attribute names extracted)
    """
    try:
        from shapely.geometry import shape
        from shapely import STRtree
    except ImportError:
        raise ImportError(
            'shapely is required for tectonic attribute assignment.  '
            'Install it with:  pip install shapely'
        ) from None

    sf = _open_shapefile(dataset)
    field_names = [f[0] for f in sf.fields[1:]]

    if dataset == 'plates':
        keep = _PLATE_FIELDS
    elif dataset in ('provinces', 'global_gprv'):
        keep = _PROVINCE_FIELDS
    else:
        keep = field_names

    geoms = []
    records = []
    for rec, sh in zip(sf.iterRecords(), sf.iterShapes()):
        if sh.shapeType == 0:
            continue
        try:
            g = shape(sh.__geo_interface__)
        except Exception:
            continue
        if g.is_empty:
            continue
        geoms.append(g)
        records.append({k: rec[field_names.index(k)]
                        for k in keep if k in field_names})

    tree = STRtree(geoms)
    return tree, geoms, records, keep


def assign_tectonic_attributes(
    lon,
    lat,
    datasets: list[str] | None = None,
) -> 'pd.DataFrame':
    """Assign tectonic attributes to geographic point locations.

    Uses point-in-polygon testing (shapely STRtree) against the tectonic
    shapefiles.  Points that fall outside all polygons receive empty strings
    for string fields and NaN for numeric fields.

    Parameters
    ----------
    lon, lat : array-like
        Decimal-degree coordinates.  NaN values are skipped and produce
        empty/NaN output rows.
    datasets : list of str, optional
        Which datasets to query.  Default ``['plates', 'provinces']``.
        Valid names: ``'plates'``, ``'provinces'``.

    Returns
    -------
    pandas.DataFrame
        Same number of rows as the input, with columns:

        From ``'plates'``: ``plate``, ``plate_code``, ``subplate``,
        ``plate_type``, ``crust_type``

        From ``'provinces'``: ``prov_name``, ``prov_type``, ``prov_group``,
        ``lastorogen``, ``continent``

    Notes
    -----
    The spatial index is built once per dataset and cached for the lifetime
    of the Python session.  Expect ~5–30 s for a first call on 1 M points.
    """
    import pandas as pd
    try:
        import shapely
    except ImportError:
        raise ImportError(
            'shapely is required for tectonic attribute assignment.'
        ) from None

    if datasets is None:
        datasets = ['plates', 'provinces']

    lon = np.asarray(lon, dtype=float).ravel()
    lat = np.asarray(lat, dtype=float).ravel()
    n   = len(lon)

    result = pd.DataFrame(index=range(n))

    for dataset in datasets:
        tree, geoms, records, fields = _build_spatial_index(dataset)

        # Default empty values
        defaults = {f: '' for f in fields}
        col_data = {f: [''] * n for f in fields}

        valid = np.where(np.isfinite(lon) & np.isfinite(lat))[0]
        if len(valid) == 0:
            result = pd.concat([result, pd.DataFrame(col_data)], axis=1)
            continue

        # Build point array (shapely 2.x vectorised constructor)
        pts = shapely.points(lon[valid], lat[valid])

        # STRtree.query returns (input_geometry_indices, tree_geometry_indices).
        # In shapely 2.x predicate='within' tests input_geom.within(tree_geom),
        # i.e. point is within polygon — which is the correct point-in-polygon test.
        pt_idx, poly_idx = tree.query(pts, predicate='within')

        # pt_idx  : indices into `pts` (0 … len(valid)-1)
        # poly_idx: indices into `geoms` / `records`
        # valid[pt_idx] maps back to original row indices
        assigned = {}   # original_row_index → record dict
        for qi, pi in zip(pt_idx, poly_idx):
            orig = int(valid[qi])
            if orig not in assigned:
                assigned[orig] = records[pi]

        for row_idx, rec in assigned.items():
            for f in fields:
                col_data[f][row_idx] = rec.get(f, '')

        result = pd.concat([result, pd.DataFrame(col_data)], axis=1)

    return result.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Polygon colour palette
# ---------------------------------------------------------------------------

def _pastel_palette(n: int) -> list[tuple]:
    """Generate *n* visually distinct pastel RGBA colours."""
    colors = []
    for i in range(n):
        h = i / n
        r, g, b = colorsys.hls_to_rgb(h, 0.88, 0.55)
        colors.append((r, g, b, 1.0))
    return colors


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def plotboundaries(
    name_or_path: str | pathlib.Path = 'plate_boundaries',
    projection: str = 'geographic',
    lon0: float = 0.0,
    lat0: float = 0.0,
    attribute: str | None = None,
    value=None,
    filters: dict | None = None,
    ax: plt.Axes | None = None,
    # barb options
    barbs: bool = False,
    barb_spacing: float | None = None,
    barb_length: float | None = None,
    barb_side: int = -1,
    barb_color=None,
    **kw,
) -> plt.Axes:
    """Plot tectonic or geographic boundaries from a shapefile.

    Parameters
    ----------
    name_or_path : str or Path
        Registered dataset name or path to a ``.shp`` file.
        Default ``'plate_boundaries'``.
    projection : str
        Projection name (same as :func:`~coast.plotcoast`).
    lon0, lat0 : float
        Central meridian and parallel in degrees.
    attribute, value : optional
        Filter to features where ``attribute == value``.
    ax : matplotlib.axes.Axes, optional
        Target axes.  Created if *None*.
    barbs : bool
        Draw filled triangular barbs perpendicular to the line (default
        ``False``).  Suitable for subduction zones.
    barb_spacing : float, optional
        Distance between barbs in map units.  Auto-scaled if *None*.
    barb_length : float, optional
        Barb length in map units.  Auto-scaled if *None*.
    barb_side : int
        ``-1`` (default) — barbs point to the left of the line direction.
        ``1`` — barbs point to the right.
    barb_color : color, optional
        Barb fill colour.  Defaults to the same colour as the line.
    **kw
        Forwarded to :func:`matplotlib.axes.Axes.plot`.

    Returns
    -------
    matplotlib.axes.Axes
    """
    kw.setdefault('color', 'C1')
    kw.setdefault('linewidth', 0.5)

    lon, lat = load_shapes(name_or_path, attribute=attribute, value=value,
                           filters=filters)
    if lon.size == 0:
        warnings.warn(f"No features loaded from '{name_or_path}'.")
        if ax is None:
            _, ax = plt.subplots()
        return ax

    lon, lat = _preprocess_lines(lon, lat, projection, lon0)
    x, y = _project(lon, lat, projection, lon0, lat0)
    x, y = _clip_lines(x, y, projection)

    if ax is None:
        _, ax = plt.subplots()

    ax.plot(x, y, **kw)

    if barbs:
        bc = barb_color if barb_color is not None else kw['color']
        _draw_barbs(ax, x, y,
                    spacing=barb_spacing,
                    length=barb_length,
                    side=barb_side,
                    color=bc)

    return ax


def plotpolygons(
    name_or_path: str | pathlib.Path = 'plates',
    projection: str = 'geographic',
    lon0: float = 0.0,
    lat0: float = 0.0,
    attribute: str | None = None,
    value=None,
    group_by: str | None = None,
    colors=None,
    alpha: float = 0.35,
    ax: plt.Axes | None = None,
    zorder: int = 0,
    **kw,
) -> plt.Axes:
    """Fill tectonic plate or province polygons on a map.

    Polygons are filled in projected space.  They are intended to be drawn
    **before** the coastline (lower *zorder*) so the coastline remains
    visible on top.

    Parameters
    ----------
    name_or_path : str or Path
        Registered dataset name or path to a polygon ``.shp`` file.
        Default ``'plates'``.
    projection : str
        Projection name matching the underlying map.
    lon0, lat0 : float
        Central meridian and parallel (must match the base map).
    attribute, value : optional
        Filter to features where ``attribute == value``.
    group_by : str, optional
        Shapefile field name to colour by.  All polygons belonging to the
        same group receive the same pastel colour.  If *None*, colours
        cycle over each individual polygon part.
    colors : list or None
        Explicit colour list.  If *None*, a pastel palette is generated
        automatically.
    alpha : float
        Fill transparency.  Default 0.35.
    ax : matplotlib.axes.Axes, optional
        Target axes.  Created if *None*.
    zorder : int
        Drawing order.  Default 0 (below coastline and boundaries).
    **kw
        Additional keyword arguments forwarded to
        :func:`matplotlib.axes.Axes.fill` (e.g. ``edgecolor``).

    Returns
    -------
    matplotlib.axes.Axes
    """
    kw.setdefault('edgecolor', 'none')
    kw.setdefault('linewidth', 0)

    if ax is None:
        _, ax = plt.subplots()

    parts = load_polygon_parts(
        name_or_path, attribute=attribute, value=value, group_by=group_by
    )

    if group_by is not None:
        # parts is a dict: group_key → list of (lon, lat) arrays
        groups = list(parts.keys())
        if colors is None:
            # Use QGIS colours for the 'plate' field; fall back to pastel palette
            if group_by == 'plate':
                palette = [
                    _QGIS_PLATE_COLORS.get(g, _pastel_palette(1)[0])
                    for g in groups
                ]
            else:
                palette = _pastel_palette(len(groups))
        else:
            palette = list(colors)
        color_map = {g: palette[i % len(palette)] for i, g in enumerate(groups)}

        for group_key, ring_list in parts.items():
            c = color_map[group_key]
            _fill_rings(ax, ring_list, projection, lon0, lat0, c, alpha, zorder, kw)
    else:
        # parts is a flat list of (lon, lat) arrays — cycle colours
        n = max(len(parts), 1)
        palette = _pastel_palette(n) if colors is None else list(colors)

        for i, (lon_ring, lat_ring) in enumerate(parts):
            c = palette[i % len(palette)]
            _fill_rings(ax, [(lon_ring, lat_ring)], projection, lon0, lat0,
                        c, alpha, zorder, kw)

    return ax


_WATERMAN_FOLDS = [-180, -90, 0, 90, 180]
_ANTIMERIDIAN    = [-180, 180]


def _densify_ring(lon, lat, max_deg=2.0):
    """Insert intermediate vertices along segments longer than max_deg degrees.

    Curved projections (e.g. Robinson) map the antimeridian to a curved edge,
    so a polygon segment with only two points at ±180° gets drawn as a straight
    line in projected space.  Densifying ensures the filled boundary follows
    the true projection curve.
    """
    out_lon = [lon[0]]
    out_lat = [lat[0]]
    for i in range(len(lon) - 1):
        dist = max(abs(lon[i + 1] - lon[i]), abs(lat[i + 1] - lat[i]))
        if dist > max_deg:
            n = int(np.ceil(dist / max_deg))
            ts = np.linspace(0, 1, n + 1)[1:]
            out_lon.extend(lon[i] + ts * (lon[i + 1] - lon[i]))
            out_lat.extend(lat[i] + ts * (lat[i + 1] - lat[i]))
        else:
            out_lon.append(lon[i + 1])
            out_lat.append(lat[i + 1])
    return np.array(out_lon), np.array(out_lat)



def _fill_rings(ax, ring_list, projection, lon0, lat0, color, alpha, zorder, kw):
    """Project and fill polygon rings.

    For Waterman projections, uses shapely polygon intersection to cut each ring
    into the 8 geographic octant panels (4 longitude bands x 2 hemispheres).
    Fold meridian splitting (via S-H) is done first to avoid shapely issues with
    antimeridian-crossing rings.  Shapely intersection then correctly inserts all
    panel boundary coordinates including panel corner vertices (e.g. the V-notch
    tip at (fold_lon, lat=0)).  Southern panels use lat=-EQ_RETRACT instead of
    lat=0 as the equatorial boundary so boundary vertices project to the southern
    Waterman face.  Fold boundary vertices are retracted slightly inside each band
    to avoid the Waterman fold-seam projection ambiguity.  The equatorial boundary
    edge is densified so fills trace the curved wing arc rather than a chord.

    For Cassini, uses the same shapely-intersection approach with 3 panels in the
    pre-rotated lon_rot frame: left back-face [-180, -90], front face [-90, +90],
    right back-face [+90, +180].  No equatorial split is needed.  The fold
    meridians at lon_rot = ±90° are curved arcs in Cassini space so boundary
    edges are densified.  Polar rings crossing the antimeridian fall back to S-H;
    non-polar rings use the two-copy approach.

    For all other projections, the ring is pre-processed with _preprocess_lines
    (NaN-based antimeridian/fold splitting) and each NaN-separated segment is
    filled independently.
    """
    key = projection.strip().lower()
    is_waterman = key in ('waterman', 'waterman_m')
    is_cassini  = key == 'cassini'

    if is_waterman:
        _FOLD_RETRACT = 1e-4
        _EQ_RETRACT   = 1e-4

        # Always clip in the unrotated (lon_rot = lon - lon0) frame so that
        # panels are always [-180,-90], [-90,0], [0,90], [90,180] regardless
        # of lon0.  This avoids (a) the sign error in computing shifted fold
        # positions, and (b) the wrapping outer band that arises when lon0≠0.
        # After clipping we project with lon0=0 since the rotation has already
        # been baked into the coordinates.
        inner_folds = [-90.0, 0.0, 90.0]
        lon_bounds  = [-180.0, -90.0, 0.0, 90.0, 180.0]

    for lon_ring, lat_ring in ring_list:
        if is_waterman:
            # Pre-rotate into the map frame: lon_rot = lon - lon0 (mod ±180).
            # Antimeridian crossings are then detected in this rotated space
            # so the outer edge of the butterfly (always at lon_rot = ±180)
            # is treated correctly for any lon0.
            lon_rot = (
                (np.asarray(lon_ring, dtype=float) - lon0 + 180.0) % 360.0
            ) - 180.0
            lat_arr = np.asarray(lat_ring, dtype=float)

            dlon = np.diff(lon_rot)
            if np.any(np.abs(dlon) > 180.0):
                # Ring crosses the rotated antimeridian (lon_rot = ±180).
                #
                # Polar rings (those that encircle a pole) must use S-H
                # pre-splitting: the two-copy approach splits them at lon=0
                # and produces geometrically inverted half-polygons that
                # fail to include the polar interior.
                #
                # Non-polar rings (simple antimeridian crossings) use the
                # two-copy approach:
                #   left  (bands 0-1): shift values > 0 by −360 so the
                #          closing edge crosses lon_rot = −180 at the outer
                #          boundary, giving shapely the correct panel edge.
                #   right (bands 2-3): shift values < 0 by +360 symmetrically.
                is_polar = lat_arr.min() < -70.0 or lat_arr.max() > 70.0
                if is_polar:
                    input_rings = [
                        (b_lon, b_lat, (0, 1, 2, 3))
                        for b_lon, b_lat
                        in _clip_ring_at_folds(lon_rot, lat_arr, inner_folds)
                    ]
                else:
                    lon_left  = lon_rot - 360.0 * (lon_rot > 0.0)
                    lon_right = lon_rot + 360.0 * (lon_rot < 0.0)
                    input_rings = [
                        (lon_left,  lat_arr, (0, 1)),
                        (lon_right, lat_arr, (2, 3)),
                    ]
            else:
                input_rings = [(lon_rot, lat_arr, (0, 1, 2, 3))]

            for b_lon, b_lat, band_ids in input_rings:
                if len(b_lon) < 3:
                    continue

                # Build shapely polygon in pre-rotated lon_rot space.
                try:
                    ring_poly = _ShpPolygon(zip(b_lon, b_lat))
                    if not ring_poly.is_valid:
                        ring_poly = ring_poly.buffer(0)
                    if ring_poly.is_empty:
                        continue
                except Exception:
                    continue

                # Intersect with the relevant panel boxes.
                for band_idx in band_ids:
                    lon_min = lon_bounds[band_idx]
                    lon_max = lon_bounds[band_idx + 1]
                    for lat_min, lat_max in ((-90.0, 0.0), (0.0, 90.0)):
                        # Southern equatorial boundary at -EQ_RETRACT so those
                        # vertices project onto the southern Waterman face.
                        clip_lat_max = lat_max if lat_max > 0.0 else -_EQ_RETRACT
                        panel_box = _shp_box(lon_min, lat_min, lon_max, clip_lat_max)

                        try:
                            inter = ring_poly.intersection(panel_box)
                        except Exception:
                            continue

                        if inter.is_empty:
                            continue

                        # Collect all Polygon pieces from the result.
                        if inter.geom_type == 'Polygon':
                            geoms = [inter]
                        elif inter.geom_type in ('MultiPolygon', 'GeometryCollection'):
                            geoms = [g for g in inter.geoms
                                     if g.geom_type == 'Polygon' and not g.is_empty]
                        else:
                            continue

                        for geom in geoms:
                            coords = np.array(geom.exterior.coords)
                            if len(coords) < 3:
                                continue
                            ilon = coords[:, 0].copy()
                            ilat = coords[:, 1]

                            # Retract fold-boundary vertices so the Waterman
                            # projection assigns them to the correct face.
                            ilon[ilon <= lon_min] = lon_min + _FOLD_RETRACT
                            ilon[ilon >= lon_max] = lon_max - _FOLD_RETRACT

                            # Densify all panel boundary edges that are curved
                            # in Waterman space: equatorial (horizontal) and
                            # fold/antimeridian (vertical).  Without densification
                            # any boundary segment projects as a straight chord
                            # rather than the actual curved arc.
                            # Northern panels: bottom edge is at lat=0 (curved).
                            # Southern panels: top edge is at clip_lat_max=-EQ_RETRACT.
                            h_dens = []
                            if lat_min == 0.0:
                                h_dens.append(0.0)
                            if clip_lat_max < 90.0:
                                h_dens.append(clip_lat_max)
                            v_dens = [lon_min + _FOLD_RETRACT, lon_max - _FOLD_RETRACT]
                            ilon, ilat = _densify_panel_boundary(
                                ilon, ilat,
                                v_bounds=v_dens, h_bounds=h_dens, max_deg=0.5,
                            )

                            # Coordinates are already in lon_rot space; project
                            # with lon0=0 (rotation was pre-applied above).
                            x, y = _project(ilon, ilat, projection, 0.0, lat0)
                            fin = np.isfinite(x) & np.isfinite(y)
                            if fin.sum() >= 3:
                                ax.fill(x[fin], y[fin], color=color, alpha=alpha,
                                        zorder=zorder, **kw)
        elif is_cassini:
            # ── Cassini panel clip ─────────────────────────────────────────
            # Pre-rotate into lon_rot = lon - lon0 frame.  Three panels:
            #   left back-face  [-180, -90]  (far hemisphere, left half)
            #   front face      [ -90,  +90] (near hemisphere)
            #   right back-face [ +90, +180] (far hemisphere, right half)
            # Fold meridians at lon_rot = ±90° are curved in Cassini space,
            # so boundary edges are densified.  Project with lon0=0 since
            # rotation is already baked in.
            _FOLD_RETRACT = 1e-4
            lon_rot = (
                (np.asarray(lon_ring, dtype=float) - lon0 + 180.0) % 360.0
            ) - 180.0
            lat_arr = np.asarray(lat_ring, dtype=float)

            dlon = np.diff(lon_rot)
            if np.any(np.abs(dlon) > 180.0):
                # Ring crosses the rotated antimeridian (lon_rot = ±180°).
                # S-H cannot be used here: when applied to find the ±90° fold
                # crossing it interpolates across the ±180° jump, creating
                # false boundary vertices for any ring — including polar ones.
                # Instead, create two continuous copies that eliminate the jump:
                #   left  copy: positive lon_rot values shifted by −360
                #   right copy: negative lon_rot values shifted by +360
                # Both copies are intersected with all three panels; panels
                # whose lon range doesn't overlap the copy simply return empty.
                # Shapely correctly handles polar interiors (e.g. Antarctica)
                # because each panel is an unambiguous lat/lon rectangle and
                # the south/north pole is naturally inside the clipped polygon.
                lon_left  = lon_rot - 360.0 * (lon_rot > 0.0)
                lon_right = lon_rot + 360.0 * (lon_rot < 0.0)
                input_rings = [
                    (lon_left,  lat_arr, (0, 1, 2)),
                    (lon_right, lat_arr, (0, 1, 2)),
                ]
            else:
                input_rings = [(lon_rot, lat_arr, (0, 1, 2))]

            # Panel boxes: idx 0 = left back, 1 = front, 2 = right back
            _cass_panels = [(-180.0, -90.0), (-90.0, 90.0), (90.0, 180.0)]

            for b_lon, b_lat, panel_ids in input_rings:
                if len(b_lon) < 3:
                    continue
                try:
                    ring_poly = _ShpPolygon(zip(b_lon, b_lat))
                    if not ring_poly.is_valid:
                        ring_poly = ring_poly.buffer(0)
                    if ring_poly.is_empty:
                        continue
                except Exception:
                    continue

                for pid in panel_ids:
                    lon_min, lon_max = _cass_panels[pid]
                    # Back-face panels (pid 0 and 2): avoid placing a vertex
                    # at exactly lat=-90.  arctan2(-inf, cos(Δlon)) = -π/2
                    # regardless of Δlon, so the south pole always projects to
                    # y=-Rπ (bottom) in the Cassini formula — but the limit as
                    # lat→-90 from above correctly wraps to y→+Rπ (top of the
                    # back face).  Keeping lat slightly above -90 sidesteps the
                    # discontinuity and lets the fill reach the top of the strip.
                    lat_min_panel = (-90.0 + _FOLD_RETRACT) if pid != 1 else -90.0
                    panel_box = _shp_box(lon_min, lat_min_panel, lon_max, 90.0)
                    try:
                        inter = ring_poly.intersection(panel_box)
                    except Exception:
                        continue
                    if inter.is_empty:
                        continue

                    if inter.geom_type == 'Polygon':
                        geoms = [inter]
                    elif inter.geom_type in ('MultiPolygon', 'GeometryCollection'):
                        geoms = [g for g in inter.geoms
                                 if g.geom_type == 'Polygon' and not g.is_empty]
                    else:
                        continue

                    for geom in geoms:
                        coords = np.array(geom.exterior.coords)
                        if len(coords) < 3:
                            continue
                        ilon = coords[:, 0].copy()
                        ilat = coords[:, 1]

                        # Retract vertices exactly on the fold seam so the
                        # Cassini projection assigns them to the correct face.
                        ilon[ilon <= lon_min] = lon_min + _FOLD_RETRACT
                        ilon[ilon >= lon_max] = lon_max - _FOLD_RETRACT

                        # Densify fold boundaries (curved arcs in Cassini space).
                        v_dens = [lon_min + _FOLD_RETRACT, lon_max - _FOLD_RETRACT]
                        ilon, ilat = _densify_panel_boundary(
                            ilon, ilat,
                            v_bounds=v_dens, h_bounds=[], max_deg=0.5,
                        )

                        x, y = _project(ilon, ilat, projection, 0.0, lat0)
                        fin = np.isfinite(x) & np.isfinite(y)
                        if fin.sum() >= 3:
                            ax.fill(x[fin], y[fin], color=color, alpha=alpha,
                                    zorder=zorder, **kw)

        else:
            lon_ring, lat_ring = _preprocess_lines(lon_ring, lat_ring, projection, lon0)
            x, y = _project(lon_ring, lat_ring, projection, lon0, lat0)

            nan_mask = ~np.isfinite(x) | ~np.isfinite(y)
            split_at = np.where(nan_mask)[0]
            starts   = np.concatenate([[0],      split_at + 1])
            ends     = np.concatenate([split_at, [len(x)]])
            for s, e in zip(starts, ends):
                seg_x = x[s:e]
                seg_y = y[s:e]
                fin   = np.isfinite(seg_x) & np.isfinite(seg_y)
                if fin.sum() >= 3:
                    ax.fill(seg_x[fin], seg_y[fin], color=color, alpha=alpha,
                            zorder=zorder, **kw)
