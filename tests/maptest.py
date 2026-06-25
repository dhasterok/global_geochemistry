"""
Smoke tests for src/maptools — visual inspection of all supported projections.

Run interactively (shows figures):
    python tests/maptest.py

Run non-interactively (saves PNGs to tests/output/, no display):
    python tests/maptest.py --save

Each figure is saved to tests/output/ regardless of mode; the directory is
created if it does not exist.  A summary line is printed for every figure so
any exception is immediately visible.
"""

import sys
import pathlib
import argparse

_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT))

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

parser = argparse.ArgumentParser(description='maptools smoke tests')
parser.add_argument('--save', action='store_true',
                    help='save PNGs and exit without displaying')
args = parser.parse_args()

if args.save:
    matplotlib.use('Agg')

from src.maptools.coast import plotcoast, tissot
from src.maptools.tectonic import plotboundaries, plotpolygons

_OUT = pathlib.Path(__file__).parent / 'output'
_OUT.mkdir(exist_ok=True)

_errors = []


def _save(fig, name):
    path = _OUT / f'{name}.png'
    fig.savefig(path, dpi=120, bbox_inches='tight')
    print(f'  ok  {name}.png')
    plt.close(fig)


def _run(label, fn):
    print(f'[ {label} ]')
    try:
        fig = fn()
        _save(fig, label)
    except Exception as exc:
        print(f'  FAIL  {exc}')
        _errors.append((label, exc))


# ── 1. Geographic ─────────────────────────────────────────────────────────────
def fig_geographic():
    fig, ax = plt.subplots(figsize=(10, 5))
    plotcoast(ax=ax)
    tissot(ax=ax)
    ax.set_title('Geographic  (lon₀ = 0°)')
    fig.tight_layout()
    return fig

_run('01_geographic', fig_geographic)


# ── 2. Geographic — Pacific centred ───────────────────────────────────────────
def fig_geographic_pacific():
    fig, ax = plt.subplots(figsize=(10, 5))
    plotcoast(lon0=180, ax=ax)
    tissot(lon0=180, ax=ax)
    ax.set_title('Geographic  (lon₀ = 180° — Pacific centred)')
    fig.tight_layout()
    return fig

_run('02_geographic_pacific', fig_geographic_pacific)


# ── 3a. Cassini  standard ─────────────────────────────────────────────────────
def fig_cassini():
    fig, ax = plt.subplots(figsize=(4, 8))
    plotcoast(projection='cassini', lon0=-70, ax=ax)
    tissot(projection='cassini', lon0=-70, ax=ax)
    ax.set_title('Cassini  (lon₀ = −70°,  lat₀ = 0°)')
    fig.tight_layout()
    return fig

_run('03a_cassini', fig_cassini)


# ── 3b. Cassini  lat0=90 ──────────────────────────────────────────────────────
def fig_cassini_lat90():
    fig, ax = plt.subplots(figsize=(4, 8))
    plotcoast(projection='cassini', lon0=-70, lat0=90, ax=ax)
    tissot(projection='cassini', lon0=-70, lat0=90, ax=ax)
    ax.set_title('Cassini  (lon₀ = −70°,  lat₀ = 90° — S pole top & bottom)')
    fig.tight_layout()
    return fig

_run('03b_cassini_lat90', fig_cassini_lat90)


# ── 4. Spilhaus ───────────────────────────────────────────────────────────────
def fig_spilhaus():
    fig, ax = plt.subplots(figsize=(7, 7))
    plotcoast(projection='spilhaus', ax=ax)
    tissot(projection='spilhaus', ax=ax)
    ax.set_title('Spilhaus  (ocean-centric)')
    fig.tight_layout()
    return fig

_run('04_spilhaus', fig_spilhaus)


# ── 5. Waterman Butterfly ─────────────────────────────────────────────────────
def fig_waterman():
    fig, ax = plt.subplots(figsize=(12, 6))
    plotcoast(projection='waterman', ax=ax)
    tissot(projection='waterman', ax=ax)
    ax.set_title('Waterman Butterfly')
    fig.tight_layout()
    return fig

_run('05_waterman_butterfly', fig_waterman)


# ── 6. Waterman M ─────────────────────────────────────────────────────────────
def fig_waterman_m():
    fig, ax = plt.subplots(figsize=(12, 5))
    plotcoast(projection='waterman_m', ax=ax)
    tissot(projection='waterman_m', ax=ax)
    ax.set_title('Waterman M')
    fig.tight_layout()
    return fig

_run('06_waterman_m', fig_waterman_m)


# ── 7. Mollweide ─────────────────────────────────────────────────────────────
def fig_mollweide():
    fig, ax = plt.subplots(figsize=(10, 5))
    plotcoast(projection='mollweide', ax=ax)
    tissot(projection='mollweide', ax=ax)
    ax.set_title("Mollweide")
    fig.tight_layout()
    return fig

_run('07_mollweide', fig_mollweide)


# ── 8. Robinson ───────────────────────────────────────────────────────────────
def fig_robinson():
    fig, ax = plt.subplots(figsize=(10, 5))
    plotcoast(projection='robinson', ax=ax)
    tissot(projection='robinson', ax=ax)
    ax.set_title("Robinson")
    fig.tight_layout()
    return fig

_run('08_robinson', fig_robinson)


# ── 9. LCC — North America ────────────────────────────────────────────────────
def fig_lcc_na():
    fig, ax = plt.subplots(figsize=(8, 7))
    plotcoast(projection='lcc', lon0=-100, lat0=40, ax=ax)
    tissot(projection='lcc', lon0=-100, lat0=40, ax=ax)
    ax.set_xlim(-5e6, 5e6)
    ax.set_ylim(-3e6, 5e6)
    ax.set_title('Lambert Conformal Conic  (lon₀ = −100°, lat₀ = 40°)')
    fig.tight_layout()
    return fig

_run('09_lcc_north_america', fig_lcc_na)


# ── 10. LCC SH — South America ────────────────────────────────────────────────
def fig_lcc_sa():
    fig, ax = plt.subplots(figsize=(7, 8))
    plotcoast(projection='lcc_sh', lon0=-60, lat0=-30, ax=ax)
    tissot(projection='lcc_sh', lon0=-60, lat0=-30, ax=ax)
    ax.set_xlim(-4e6, 4e6)
    ax.set_ylim(-5e6, 3e6)
    ax.set_title('LCC Southern Hemisphere  (lon₀ = −60°)')
    fig.tight_layout()
    return fig

_run('10_lcc_south_america', fig_lcc_sa)


# ── 11. Albers EA — Europe ────────────────────────────────────────────────────
def fig_aea_europe():
    fig, ax = plt.subplots(figsize=(9, 7))
    plotcoast(projection='aea', lon0=10, lat0=50, ax=ax)
    tissot(projection='aea', lon0=10, lat0=50, ax=ax)
    ax.set_xlim(-4e6, 4e6)
    ax.set_ylim(-3e6, 4e6)
    ax.set_title('Albers Equal-Area Conic  (lon₀ = 10°, lat₀ = 50°)')
    fig.tight_layout()
    return fig

_run('11_aea_europe', fig_aea_europe)


# ── 12. Albers EA SH — Australia ──────────────────────────────────────────────
def fig_aea_au():
    fig, ax = plt.subplots(figsize=(9, 7))
    plotcoast(projection='aea_sh', lon0=135, lat0=-30, ax=ax)
    tissot(projection='aea_sh', lon0=135, lat0=-30, ax=ax)
    ax.set_xlim(-4e6, 4e6)
    ax.set_ylim(-4e6, 3e6)
    ax.set_title('Albers EA Southern Hemisphere  (lon₀ = 135°)')
    fig.tight_layout()
    return fig

_run('12_aea_australia', fig_aea_au)


# ── 13. Tectonic boundaries — Robinson global ─────────────────────────────────
def fig_tectonic_global():
    fig, ax = plt.subplots(figsize=(12, 6))
    plotcoast(projection='robinson', ax=ax)
    plotboundaries('plate_boundaries', projection='robinson', ax=ax,
                   attribute='type', value='spreading center',
                   color='C3', linewidth=0.8)
    plotboundaries('plate_boundaries', projection='robinson', ax=ax,
                   attribute='type', value='subduction zone',
                   color='C0', linewidth=0.8)
    ax.set_title('Robinson — plate boundaries  (red: spreading, blue: subduction)')
    ax.legend(handles=[
        mpatches.Patch(color='C3', label='spreading center'),
        mpatches.Patch(color='C0', label='subduction zone'),
    ], fontsize=8, loc='lower left')
    fig.tight_layout()
    return fig

_run('13_tectonic_robinson', fig_tectonic_global)


# ── 14. Tectonic boundaries — LCC regional ────────────────────────────────────
def fig_tectonic_lcc():
    fig, ax = plt.subplots(figsize=(9, 7))
    plotcoast(projection='lcc', lon0=-100, lat0=40, ax=ax)
    plotboundaries('plate_boundaries', projection='lcc', lon0=-100, lat0=40,
                   ax=ax, color='C3', linewidth=0.8)
    plotboundaries('plates', projection='lcc', lon0=-100, lat0=40,
                   ax=ax, color='C0', linewidth=0.4, alpha=0.4)
    ax.set_xlim(-5e6, 5e6)
    ax.set_ylim(-3e6, 5e6)
    ax.set_title('LCC — North America  (plate boundaries + plate polygons)')
    fig.tight_layout()
    return fig

_run('14_tectonic_lcc', fig_tectonic_lcc)


# ── 15. Waterman — tectonic overlay ───────────────────────────────────────────
def fig_tectonic_waterman():
    fig, ax = plt.subplots(figsize=(13, 7))
    plotcoast(projection='waterman', ax=ax)
    plotboundaries('plate_boundaries', projection='waterman', ax=ax,
                   color='C3', linewidth=0.7)
    ax.set_title('Waterman Butterfly — plate boundaries')
    fig.tight_layout()
    return fig

_run('15_tectonic_waterman', fig_tectonic_waterman)


# ── 16. Plate polygons — Robinson, pastel fill below coastline ────────────────
def fig_plate_polygons_robinson():
    fig, ax = plt.subplots(figsize=(12, 6))
    plotpolygons('plates', projection='robinson', ax=ax, zorder=0)
    plotcoast(projection='robinson', ax=ax)
    plotboundaries('plate_boundaries', projection='robinson', ax=ax,
                   color='0.3', linewidth=0.5)
    ax.set_title('Robinson — plate polygons (pastel) + boundaries')
    fig.tight_layout()
    return fig

_run('16_plate_polygons_robinson', fig_plate_polygons_robinson)


# ── 17. Plate polygons — geographic, grouped by plate ────────────────────────
def fig_plate_polygons_grouped():
    fig, ax = plt.subplots(figsize=(12, 5))
    plotpolygons('plates', group_by='plate', ax=ax, zorder=0)
    plotcoast(ax=ax)
    plotboundaries('plate_boundaries', ax=ax, color='0.4', linewidth=0.4)
    ax.set_title('Geographic — plate polygons coloured by plate name')
    fig.tight_layout()
    return fig

_run('17_plate_polygons_grouped', fig_plate_polygons_grouped)


# ── 18. Subduction zone barbs — Robinson global ───────────────────────────────
def fig_barbs_robinson():
    fig, ax = plt.subplots(figsize=(12, 6))
    plotpolygons('plates', projection='robinson', ax=ax, zorder=0, alpha=0.25)
    plotcoast(projection='robinson', ax=ax)
    plotboundaries('plate_boundaries', projection='robinson', ax=ax,
                   attribute='type', value='spreading center',
                   color='C3', linewidth=0.8)
    plotboundaries('plate_boundaries', projection='robinson', ax=ax,
                   attribute='type', value='subduction zone',
                   color='k', linewidth=0.8,
                   barbs=True, barb_side=-1)
    ax.set_title('Robinson — subduction barbs (side=1) + spreading centres')
    ax.legend(handles=[
        mpatches.Patch(color='C3', label='spreading center'),
        mpatches.Patch(color='k',  label='subduction zone'),
    ], fontsize=8, loc='lower left')
    fig.tight_layout()
    return fig

_run('18_barbs_robinson', fig_barbs_robinson)


# ── 19. Subduction barbs — LCC North America ─────────────────────────────────
def fig_barbs_lcc():
    fig, ax = plt.subplots(figsize=(9, 7))
    plotpolygons('plates', projection='lcc', lon0=-100, lat0=40,
                 ax=ax, zorder=0, alpha=0.3)
    plotcoast(projection='lcc', lon0=-100, lat0=40, ax=ax)
    plotboundaries('plate_boundaries', projection='lcc', lon0=-100, lat0=40,
                   ax=ax, attribute='type', value='subduction zone',
                   color='k', linewidth=1.0,
                   barbs=True, barb_side=-1)
    plotboundaries('plate_boundaries', projection='lcc', lon0=-100, lat0=40,
                   ax=ax, attribute='type', value='spreading center',
                   color='C3', linewidth=0.8)
    ax.set_xlim(-5e6, 5e6)
    ax.set_ylim(-3e6, 5e6)
    ax.set_title('LCC — North America, subduction barbs')
    fig.tight_layout()
    return fig

_run('19_barbs_lcc', fig_barbs_lcc)


# ── 20. Subduction barbs — Waterman ──────────────────────────────────────────
def fig_barbs_waterman():
    fig, ax = plt.subplots(figsize=(13, 7))
    plotpolygons('plates', projection='waterman', ax=ax, zorder=0, alpha=0.25)
    plotcoast(projection='waterman', ax=ax)
    plotboundaries('plate_boundaries', projection='waterman', ax=ax,
                   attribute='type', value='subduction zone',
                   color='k', linewidth=0.8,
                   barbs=True, barb_side=-1)
    plotboundaries('plate_boundaries', projection='waterman', ax=ax,
                   attribute='type', value='spreading center',
                   color='C3', linewidth=0.7)
    ax.set_title('Waterman — plate polygons + subduction barbs + spreading')
    fig.tight_layout()
    return fig

_run('20_barbs_waterman', fig_barbs_waterman)


# ── Summary ───────────────────────────────────────────────────────────────────
print()
if _errors:
    print(f'FAILED ({len(_errors)} error(s)):')
    for name, exc in _errors:
        print(f'  {name}: {exc}')
    sys.exit(1)
else:
    n = 20 - len(_errors)
    print(f'All {n} smoke tests passed.  PNGs saved to tests/output/')

if not args.save:
    plt.show()
