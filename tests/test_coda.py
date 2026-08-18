"""Tests for the compositional-data-analysis re-exports in geochem.coda
that LaME's HDBSCAN clustering path (LaserMapExplorer/src/app/DataAnalysis.py)
depends on: closure, multiplicative_replacement, clr, anti_clr.

Pure numpy -- no PyQt/QApplication needed.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root / "src"))

from geochem.coda import anti_clr, closure, clr, multiplicative_replacement


def test_closure_rows_sum_to_one():
    X = np.array([[1.0, 2.0, 3.0], [4.0, 0.5, 1.5]])
    closed = closure(X)
    assert np.allclose(closed.sum(axis=1), 1.0)


def test_clr_rows_sum_to_zero():
    X = closure(np.array([[1.0, 2.0, 3.0], [10.0, 20.0, 30.0], [5.0, 5.0, 5.0]]))
    Y = clr(X)
    assert np.allclose(Y.sum(axis=1), 0.0, atol=1e-10)


def test_clr_scale_invariant_per_row():
    """clr(close(c * x)) == clr(close(x)) for any positive per-row scaling c
    -- the defining property that makes clr appropriate for compositional
    data whose absolute total (e.g. total ablated mass) isn't meaningful.
    """
    X = np.array([[1.0, 4.0, 20.0], [3.0, 3.0, 3.0]])
    scaled = X * np.array([[10.0], [0.001]])
    assert np.allclose(clr(closure(X)), clr(closure(scaled)), atol=1e-10)


def test_anti_clr_round_trip():
    X = closure(np.array([[1.0, 2.0, 3.0], [10.0, 20.0, 30.0], [7.0, 1.0, 2.0]]))
    recovered = anti_clr(clr(X))
    assert np.allclose(recovered, X, atol=1e-10)


def test_clr_on_raw_zero_containing_row_is_not_finite():
    """clr alone does not handle zeros (matches composition_stats' own
    contract) -- LaME's clustering path relies on this and always runs
    ``multiplicative_replacement`` first; this test locks in that plain
    ``clr`` on an unreplaced zero is *not* silently finite/wrong.
    """
    X = closure(np.array([[1.0, 0.0, 3.0]]))
    with np.errstate(divide='ignore', invalid='ignore'):
        Y = clr(X)
    assert not np.all(np.isfinite(Y))


def test_multiplicative_replacement_then_clr_is_finite_and_preserves_other_ratios():
    X = closure(np.array([[1.0, 0.0, 3.0], [2.0, 5.0, 4.0]]))
    replaced = multiplicative_replacement(X)
    Y = clr(replaced)
    assert np.all(np.isfinite(Y))
    # the non-replaced parts' ratio in the zero-containing row is (nearly)
    # untouched by the replacement -- only the zero itself was materially changed
    assert Y[0, 0] - Y[0, 2] == pytest.approx(np.log(X[0, 0] / X[0, 2]), abs=1e-3)
