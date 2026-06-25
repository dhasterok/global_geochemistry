"""
Smoke tests for src/plotting/radar — radarplot() and radar_prep().

Run interactively (shows figures):
    python tests/radartest.py

Run non-interactively (saves PNGs to tests/output/, no display):
    python tests/radartest.py --save
"""

import sys
import pathlib
import argparse

_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT))

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='radar smoke tests')
parser.add_argument('--save', action='store_true',
                    help='save PNGs and exit without displaying')
args = parser.parse_args()

if args.save:
    matplotlib.use('Agg')

from src.plotting.radar import radar_prep, radarplot

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
        import traceback
        traceback.print_exc()
        print(f'  FAIL  {exc}')
        _errors.append((label, exc))


# ── Synthetic geochemistry data ────────────────────────────────────────────────
rng = np.random.default_rng(42)

fields = ['SiO2', 'TiO2', 'Al2O3', 'FeOT', 'MgO', 'CaO', 'Na2O', 'K2O']


def _make_df(n=120):
    import pandas as pd
    rock_params = {
        'basalt':   {'SiO2': (48, 3),  'TiO2': (1.5, 0.5), 'Al2O3': (16, 2),
                     'FeOT': (9, 2),   'MgO':  (7, 2),      'CaO':   (10, 2),
                     'Na2O': (2.5, 0.5), 'K2O': (0.5, 0.3)},
        'andesite': {'SiO2': (58, 3),  'TiO2': (0.9, 0.3), 'Al2O3': (17, 2),
                     'FeOT': (6, 2),   'MgO':  (4, 1.5),    'CaO':   (7, 2),
                     'Na2O': (3.5, 0.5), 'K2O': (1.5, 0.5)},
        'rhyolite': {'SiO2': (74, 3),  'TiO2': (0.2, 0.1), 'Al2O3': (13, 1),
                     'FeOT': (2, 1),   'MgO':  (0.5, 0.3),  'CaO':   (1.5, 0.7),
                     'Na2O': (4.0, 0.5), 'K2O': (4.5, 0.7)},
    }
    rows = []
    for rock, params in rock_params.items():
        for _ in range(n // 3):
            row = {'rock_type': rock}
            for f in fields:
                mu, sig = params[f]
                row[f] = max(0.01, rng.normal(mu, sig))
            rows.append(row)
    return pd.DataFrame(rows)


# ── 1. All axes labeled (default) ────────────────────────────────────────────
def fig_all_axes_labeled():
    df = _make_df()
    prep = radar_prep(df, fields, group_field='rock_type',
                      groups=['basalt', 'andesite', 'rhyolite'],
                      quantiles=[0.25, 0.50, 0.75])
    fig, ax = plt.subplots(figsize=(7, 7))
    radarplot(**prep, ax=ax)  # tick_fields='all' is default
    ax.set_title('Web — all axes labeled (default)', pad=20)
    fig.tight_layout()
    return fig

_run('01_all_axes_labeled', fig_all_axes_labeled)


# ── 2. No tick labels ────────────────────────────────────────────────────────
def fig_no_tick_labels():
    df = _make_df()
    prep = radar_prep(df, fields, group_field='rock_type',
                      groups=['basalt', 'andesite', 'rhyolite'],
                      quantiles=[0.25, 0.50, 0.75])
    fig, ax = plt.subplots(figsize=(7, 7))
    radarplot(**prep, ax=ax, tick_fields=None)
    ax.set_title('Web — tick_fields=None (no labels)', pad=20)
    fig.tight_layout()
    return fig

_run('02_no_tick_labels', fig_no_tick_labels)


# ── 3. Selected axes labeled ──────────────────────────────────────────────────
def fig_selected_tick_fields():
    df = _make_df()
    prep = radar_prep(df, fields, group_field='rock_type',
                      groups=['basalt', 'andesite', 'rhyolite'],
                      quantiles=[0.25, 0.50, 0.75])
    fig, ax = plt.subplots(figsize=(7, 7))
    radarplot(**prep, ax=ax, tick_fields=[0, 2, 4, 6])
    ax.set_title('Web — tick_fields=[0, 2, 4, 6] (alternating)', pad=20)
    fig.tight_layout()
    return fig

_run('03_selected_tick_fields', fig_selected_tick_fields)


# ── 4. Custom axis display labels (dict override) ─────────────────────────────
def fig_axis_labels():
    df = _make_df()
    prep = radar_prep(df, fields, group_field='rock_type',
                      groups=['basalt', 'andesite', 'rhyolite'],
                      quantiles=[0.25, 0.50, 0.75])
    fig, ax = plt.subplots(figsize=(7, 7))
    radarplot(**prep, ax=ax,
              axis_labels={'TiO2': 'TiO₂', 'Al2O3': 'Al₂O₃', 'FeOT': 'FeO*',
                           'Na2O': 'Na₂O', 'K2O': 'K₂O'})
    ax.set_title('Web — custom axis labels (dict override)', pad=20)
    fig.tight_layout()
    return fig

_run('04_custom_axis_labels', fig_axis_labels)


# ── 5. Custom axis labels (full list) ────────────────────────────────────────
def fig_axis_labels_list():
    df = _make_df()
    prep = radar_prep(df, fields, group_field='rock_type',
                      groups=['basalt', 'andesite', 'rhyolite'])
    fig, ax = plt.subplots(figsize=(7, 7))
    display = ['SiO₂', 'TiO₂', 'Al₂O₃', 'FeO*', 'MgO', 'CaO', 'Na₂O', 'K₂O']
    radarplot(**prep, ax=ax, axis_labels=display)
    ax.set_title('Web — full axis_labels list (unicode)', pad=20)
    fig.tight_layout()
    return fig

_run('05_axis_labels_list', fig_axis_labels_list)


# ── 6. Per-axis range: dict, partial override ─────────────────────────────────
def fig_axis_lim_dict():
    df = _make_df()
    prep = radar_prep(df, fields, group_field='rock_type',
                      groups=['basalt', 'andesite', 'rhyolite'],
                      quantiles=[0.25, 0.50, 0.75])
    fig, ax = plt.subplots(figsize=(7, 7))
    radarplot(**prep, ax=ax,
              axis_lim={'SiO2': (40, 80), 'MgO': (0, 15), 'CaO': (0, 15)})
    ax.set_title('Web — axis_lim dict (partial override)', pad=20)
    fig.tight_layout()
    return fig

_run('06_axis_lim_dict', fig_axis_lim_dict)


# ── 7. Per-axis range: full (2, n_fields) array ───────────────────────────────
def fig_axis_lim_array():
    df = _make_df()
    prep = radar_prep(df, fields, group_field='rock_type',
                      groups=['basalt', 'andesite', 'rhyolite'],
                      quantiles=[0.25, 0.50, 0.75])
    axis_lim = np.array([
        [40,  0, 10,  0,  0,  0,  0, 0],
        [80,  3, 20, 15, 15, 15,  6, 8],
    ])
    fig, ax = plt.subplots(figsize=(7, 7))
    radarplot(**prep, ax=ax, axis_lim=axis_lim)
    ax.set_title('Web — axis_lim (2, n_fields) array', pad=20)
    fig.tight_layout()
    return fig

_run('07_axis_lim_array', fig_axis_lim_array)


# ── 8. Log scale — selected axes by name ──────────────────────────────────────
def fig_log_named():
    df = _make_df()
    prep = radar_prep(df, fields, group_field='rock_type',
                      groups=['basalt', 'andesite', 'rhyolite'],
                      quantiles=[0.25, 0.50, 0.75])
    fig, ax = plt.subplots(figsize=(7, 7))
    radarplot(**prep, ax=ax,
              log=['TiO2', 'MgO', 'K2O'],
              axis_labels={'TiO2': 'TiO₂ (log)', 'MgO': 'MgO (log)',
                           'K2O': 'K₂O (log)'})
    ax.set_title('Web — log scale on TiO₂, MgO, K₂O', pad=20)
    fig.tight_layout()
    return fig

_run('08_log_named_axes', fig_log_named)


# ── 9. Log scale — all axes ───────────────────────────────────────────────────
def fig_log_all():
    df = _make_df()
    # Use only trace-element-style fields so log makes sense
    trace_fields = ['TiO2', 'MgO', 'K2O', 'Na2O']
    prep = radar_prep(df, trace_fields, group_field='rock_type',
                      groups=['basalt', 'andesite', 'rhyolite'],
                      quantiles=[0.25, 0.50, 0.75])
    fig, ax = plt.subplots(figsize=(7, 7))
    radarplot(**prep, ax=ax, log=True, n_ticks=3)
    ax.set_title('Web — log=True (all axes log scale)', pad=20)
    fig.tight_layout()
    return fig

_run('09_log_all_axes', fig_log_all)


# ── 10. Log scale — by index ──────────────────────────────────────────────────
def fig_log_by_index():
    df = _make_df()
    prep = radar_prep(df, fields, group_field='rock_type',
                      groups=['basalt', 'andesite', 'rhyolite'])
    fig, ax = plt.subplots(figsize=(7, 7))
    radarplot(**prep, ax=ax, log=[1, 4, 7])  # TiO2, MgO, K2O by index
    ax.set_title('Web — log=[1, 4, 7] (by index)', pad=20)
    fig.tight_layout()
    return fig

_run('10_log_by_index', fig_log_by_index)


# ── 11. Custom colours ────────────────────────────────────────────────────────
def fig_custom_colors():
    df = _make_df()
    prep = radar_prep(df, fields, group_field='rock_type',
                      groups=['basalt', 'andesite', 'rhyolite'],
                      quantiles=[0.25, 0.50, 0.75])
    fig, ax = plt.subplots(figsize=(7, 7))
    radarplot(**prep, ax=ax,
              colors=['#e6194b', '#3cb44b', '#4363d8'],
              alpha_fill=0.25)
    ax.set_title('Web — custom colors per group', pad=20)
    fig.tight_layout()
    return fig

_run('11_custom_colors', fig_custom_colors)


# ── 12. Five quantiles ────────────────────────────────────────────────────────
def fig_five_quantiles():
    df = _make_df()
    prep = radar_prep(df, fields, group_field='rock_type',
                      groups=['basalt'],
                      quantiles=[0.05, 0.25, 0.50, 0.75, 0.95])
    fig, ax = plt.subplots(figsize=(7, 7))
    radarplot(**prep, ax=ax, n_ticks=5)
    ax.set_title('Web — five quantiles (basalt)', pad=20)
    fig.tight_layout()
    return fig

_run('12_five_quantiles', fig_five_quantiles)


# ── 13. Series — clockwise ────────────────────────────────────────────────────
def fig_series_cw():
    n_obs = 60
    t = np.linspace(0, 1, n_obs)
    data = np.column_stack([
        48 + 26 * t + rng.normal(0, 1.5, n_obs),
        np.clip(1.5 - 1.3 * t + rng.normal(0, 0.2, n_obs), 0.05, None),
        16 - 3 * t + rng.normal(0, 1, n_obs),
        np.clip(9 - 7 * t + rng.normal(0, 1, n_obs), 0.1, None),
        np.clip(7 - 6.5 * t + rng.normal(0, 0.8, n_obs), 0.05, None),
        np.clip(10 - 8.5 * t + rng.normal(0, 1, n_obs), 0.05, None),
        2.5 + 1.5 * t + rng.normal(0, 0.3, n_obs),
        np.clip(0.5 + 4 * t + rng.normal(0, 0.4, n_obs), 0.05, None),
    ])
    fig, ax = plt.subplots(figsize=(7, 7))
    radarplot(vals=data, fields=fields, style='series',
              clockwise=True, cmap='plasma', ax=ax,
              log=['TiO2', 'MgO', 'CaO', 'K2O'])
    ax.set_title('Series — clockwise, log on 4 axes', pad=20)
    fig.tight_layout()
    return fig

_run('13_series_clockwise', fig_series_cw)


# ── 14. Series — counterclockwise, tick_fields=None ──────────────────────────
def fig_series_ccw():
    n_obs = 50
    t = np.linspace(0, 1, n_obs)
    data = np.column_stack([
        74 - 26 * t + rng.normal(0, 1.5, n_obs),
        np.clip(0.2 + 1.3 * t + rng.normal(0, 0.15, n_obs), 0.01, None),
        13 + 3 * t + rng.normal(0, 1, n_obs),
        np.clip(2 + 7 * t + rng.normal(0, 1, n_obs), 0.1, None),
        np.clip(0.5 + 6.5 * t + rng.normal(0, 0.5, n_obs), 0.01, None),
        np.clip(1.5 + 8.5 * t + rng.normal(0, 0.8, n_obs), 0.05, None),
        4.0 - 1.5 * t + rng.normal(0, 0.3, n_obs),
        np.clip(4.5 - 4 * t + rng.normal(0, 0.4, n_obs), 0.05, None),
    ])
    fig, ax = plt.subplots(figsize=(7, 7))
    radarplot(vals=data, fields=fields, style='series',
              clockwise=False, cmap='coolwarm', ax=ax,
              tick_fields=None)
    ax.set_title('Series — CCW, tick_fields=None', pad=20)
    fig.tight_layout()
    return fig

_run('14_series_counterclockwise', fig_series_ccw)


# ── 15. Combined: custom lim + log + axis_labels + selected ticks ─────────────
def fig_combined():
    df = _make_df()
    prep = radar_prep(df, fields, group_field='rock_type',
                      groups=['basalt', 'andesite', 'rhyolite'],
                      quantiles=[0.25, 0.50, 0.75])
    fig, ax = plt.subplots(figsize=(8, 8))
    radarplot(**prep, ax=ax,
              axis_lim={'SiO2': (40, 80), 'TiO2': (0.1, 3.0),
                        'MgO': (0.3, 15), 'K2O': (0.05, 8)},
              log=['TiO2', 'MgO', 'K2O'],
              axis_labels={'SiO2': 'SiO₂', 'TiO2': 'TiO₂ (log)',
                           'Al2O3': 'Al₂O₃', 'FeOT': 'FeO*',
                           'MgO': 'MgO (log)', 'CaO': 'CaO',
                           'Na2O': 'Na₂O', 'K2O': 'K₂O (log)'},
              tick_fields=[0, 1, 4, 7],
              colors=['steelblue', 'darkorange', 'forestgreen'],
              alpha_fill=0.2, n_ticks=4)
    ax.set_title('Combined: custom lim + log + labels + selected ticks', pad=20)
    fig.tight_layout()
    return fig

_run('15_combined', fig_combined)


# ── Summary ────────────────────────────────────────────────────────────────────
print()
if _errors:
    print(f'FAILED ({len(_errors)} error(s)):')
    for name, exc in _errors:
        print(f'  {name}: {exc}')
    sys.exit(1)
else:
    n = 15 - len(_errors)
    print(f'All {n} radar smoke tests passed.  PNGs saved to tests/output/')

if not args.save:
    plt.show()
