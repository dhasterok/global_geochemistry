"""
Custom colormaps for geochemical plotting.

Named colormaps are loaded from ``src/data/custom_colormaps.csv``.  Each
row in that file is ``name, #hex1, #hex2, ..., #hexN`` and the colours
are used as equally-spaced stops in a :class:`~matplotlib.colors.LinearSegmentedColormap`.

Three-point palatte colormaps from MATLAB ``colormap2.m`` are also
registered (``'red'``, ``'blue'``, ``'rwb'``, etc.).

Any named matplotlib colormap is accepted as a fallback.

Main API
--------
``get_colormap(name, n=256, reverse=False)``
    Return a :class:`~matplotlib.colors.LinearSegmentedColormap`.

``palatte(left, center, right, n=256)``
    Build an N×3 RGB array interpolated through three control colours.

``available_colormaps()``
    List all registered names.
"""

from __future__ import annotations

import importlib.resources
import pathlib
import re

import matplotlib.cm as mpl_cm
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, to_rgb

# ---------------------------------------------------------------------------
# Three-point colormap builder (equivalent to MATLAB palatte.m)
# ---------------------------------------------------------------------------

def palatte(left, center, right, n: int = 256) -> np.ndarray:
    """Build an N×3 RGB array interpolated left → center → right.

    Equivalent to MATLAB ``palatte()``.  The centre colour lands at the
    midpoint of the scale.  If *n* is even it is incremented by one so
    the midpoint is always an integer index.

    Parameters
    ----------
    left, center, right : array-like of length 3
        RGB triplets in [0, 1].
    n : int
        Number of colour steps (default 256).

    Returns
    -------
    numpy.ndarray, shape (n, 3)
    """
    left   = np.asarray(left,   dtype=float)
    center = np.asarray(center, dtype=float)
    right  = np.asarray(right,  dtype=float)

    if n % 2 == 0:
        n += 1

    t = np.linspace(0.0, 1.0, n)
    colors = np.where(
        t[:, None] <= 0.5,
        left   + (center - left)  * (t[:, None] * 2),
        center + (right  - center) * ((t[:, None] - 0.5) * 2),
    )
    return np.clip(colors, 0.0, 1.0)


# ---------------------------------------------------------------------------
# Registry — populated at import time
# ---------------------------------------------------------------------------

_REGISTRY: dict[str, LinearSegmentedColormap] = {}


def _register_palatte(name: str, left, center, right, n: int = 256) -> None:
    colors = palatte(left, center, right, n)
    _REGISTRY[name] = LinearSegmentedColormap.from_list(name, colors)


def _hex_to_rgb(h: str) -> tuple[float, float, float]:
    h = h.strip().lstrip('#')
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return r / 255.0, g / 255.0, b / 255.0


def _colormap_csv_path() -> pathlib.Path | None:
    """Return the path to ``custom_colormaps.csv``.

    Prefers the version shipped with the installed ``blueberry-colortools``
    package (``blueberry.COLORMAP_PATH``).  Falls back to the copy stored
    alongside this package in ``src/data/``.
    """
    try:
        import blueberry
        p = pathlib.Path(blueberry.COLORMAP_PATH) / 'custom_colormaps.csv'
        if p.exists():
            return p
    except ImportError:
        pass

    p = pathlib.Path(__file__).parent.parent / 'data' / 'custom_colormaps.csv'
    return p if p.exists() else None


def _load_csv() -> None:
    """Parse custom_colormaps.csv and register all entries."""
    csv_path = _colormap_csv_path()
    if csv_path is None:
        return
    with csv_path.open(encoding='utf-8-sig') as fh:
        for line in fh:
            parts = [p.strip() for p in line.rstrip('\n').split(',')]
            if not parts or not parts[0]:
                continue
            name = parts[0].lower()
            hexes = [p for p in parts[1:] if re.match(r'^#[0-9A-Fa-f]{6}$', p)]
            if len(hexes) < 2:
                continue
            colors = [_hex_to_rgb(h) for h in hexes]
            _REGISTRY[name] = LinearSegmentedColormap.from_list(name, colors)


# Three-point colormaps from colormap2.m
_register_palatte('red',    [0.9]*3,                [1, .5,  0],         [.8627,.0784,.2353])
_register_palatte('red2',   [1, 1, 1],              [1, .5,  0],         [.8627,.0784,.2353])
_register_palatte('blue',   [0.9]*3,                [0, .9804,.6039],    [0,  0, .8039])
_register_palatte('green',  [0.9]*3,                [1, .7843,.1176],    [.1333,.5451,.1333])
_register_palatte('violet', [0.9]*3,                [1, .6275,.4784],    [.5804,0, .8275])
_register_palatte('brown',  [0.9]*3,                [.8706,.7216,.5294], [.5451,.2706,.0745])
_register_palatte('rwb',    [.8039,  0,    0],      [1, 1, 1],           [0,  0, .8039])
_register_palatte('rwb2',   [.8039,.3608,.3608],    [1, 1, 1],           [.3922,.5843,.9294])
_register_palatte('ryb',    [.8039,  0,    0],      [1, .7843,.1176],    [0,  0, .8039])
_register_palatte('ryb2',   [.8039,.3608,.3608],    [.9412,.902,.549],   [.3922,.5843,.9294])

# Load CSV entries (may override the three-point ones above if names collide)
_load_csv()


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def get_colormap(
    name: str | np.ndarray,
    n: int = 256,
    reverse: bool = False,
) -> LinearSegmentedColormap:
    """Return a matplotlib colormap by name.

    Parameters
    ----------
    name : str or array-like, shape (3, 3)
        Named colormap or a three-row array ``[left, center, right]`` of
        RGB triplets in [0, 1] for an on-the-fly three-point palatte.
    n : int
        Number of colour levels (default 256).
    reverse : bool
        Reverse the colormap direction.

    Returns
    -------
    matplotlib.colors.LinearSegmentedColormap
    """
    if not isinstance(name, str):
        arr = np.asarray(name, dtype=float)
        if arr.shape != (3, 3):
            raise ValueError("Color control array must have shape (3, 3): [left, center, right]")
        colors = palatte(arr[0], arr[1], arr[2], n)
        cmap = LinearSegmentedColormap.from_list('custom', colors)
        return cmap.reversed() if reverse else cmap

    key = name.lower()
    if key in _REGISTRY:
        cmap = _REGISTRY[key].resampled(n)
        return cmap.reversed() if reverse else cmap

    # Fall back to cmcrameri
    try:
        import cmcrameri.cm as cmc
        for attr in (key, key.replace(' ', '_'), key.replace('-', '_')):
            if hasattr(cmc, attr):
                cmap = getattr(cmc, attr)
                return cmap.reversed() if reverse else cmap
    except ImportError:
        pass

    # Fall back to matplotlib
    try:
        cmap = mpl_cm.get_cmap(name, n)
        return cmap.reversed() if reverse else cmap
    except (ValueError, KeyError):
        raise ValueError(
            f"Unknown colormap '{name}'. "
            f"Available custom names: {available_colormaps()}.\n"
            "Also accepts any cmcrameri or matplotlib colormap name."
        )


# Convenience alias matching the shorter name used in plotting modules.
get_cmap = get_colormap


def available_colormaps() -> list[str]:
    """Return sorted list of all registered custom colormap names."""
    return sorted(_REGISTRY.keys())


def list_cmaps(source: str = 'custom') -> None:
    """Print available colormap names.

    Parameters
    ----------
    source : {'custom', 'cmcrameri', 'matplotlib', 'all'}
    """
    if source in ('custom', 'all'):
        names = available_colormaps()
        print(f"Custom ({len(names)}):")
        for n in names:
            print(f"  {n}")

    if source in ('cmcrameri', 'all'):
        try:
            import cmcrameri.cm as cmc
            names = sorted(x for x in dir(cmc)
                           if not x.startswith('_') and not x.endswith('_r'))
            print(f"\ncmcrameri ({len(names)}):")
            for n in names:
                print(f"  {n}")
        except ImportError:
            print("\ncmcrameri not installed.")

    if source in ('matplotlib', 'all'):
        import matplotlib.pyplot as plt
        names = sorted(plt.colormaps())
        print(f"\nMatplotlib ({len(names)}):")
        for n in names:
            print(f"  {n}")
