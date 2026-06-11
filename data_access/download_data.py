"""
Data access utilities for the global geochemistry project.

Google Drive automation is not yet implemented.  Use list_datasets() to
see what datasets are registered; download_dataset() raises NotImplementedError
until a download backend is wired up.
"""

from __future__ import annotations

import json
import pathlib

_REGISTRY_PATH = pathlib.Path(__file__).parent / 'data_registry.json'


def list_datasets() -> list[dict]:
    """Return all registered datasets from data_registry.json."""
    with _REGISTRY_PATH.open() as f:
        registry = json.load(f)
    return registry['datasets']


def download_dataset(name: str, dest: str | pathlib.Path | None = None) -> pathlib.Path:
    """Download a registered dataset to *dest*.

    Parameters
    ----------
    name : str
        Dataset name as listed in data_registry.json.
    dest : path-like, optional
        Destination directory.  Defaults to the repo data/ directory.

    Raises
    ------
    NotImplementedError
        Google Drive download automation is not yet implemented.
    """
    raise NotImplementedError(
        f"Automated download for '{name}' is not yet implemented. "
        "Copy the dataset manually — see data_registry.json for its location."
    )
