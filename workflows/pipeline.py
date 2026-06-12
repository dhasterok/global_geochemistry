"""
Importable pipeline entry-point for the global geochemistry database.

``run_pipeline()`` is an alias for
:func:`src.geochem.fileio.pipeline.load_geochem` kept for backward
compatibility with existing notebooks and analysis scripts.
"""

from __future__ import annotations

import sys
import pathlib

_root = pathlib.Path(__file__).resolve().parent.parent
if str(_root) not in sys.path:
    sys.path.insert(0, str(_root))

from src.geochem.fileio.pipeline import load_geochem  # noqa: E402

run_pipeline = load_geochem

__all__ = ['load_geochem', 'run_pipeline']
