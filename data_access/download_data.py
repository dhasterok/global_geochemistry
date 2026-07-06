"""
Data access utilities for the global geochemistry project.

Use :func:`list_datasets` to see registered datasets and
:func:`download_dataset` to automatically fetch those with a known URL.
After a successful download the local path is written back into
``data_registry.json`` so subsequent imports resolve without re-downloading.
"""

from __future__ import annotations

import json
import pathlib
import shutil
import urllib.request
import zipfile

_REGISTRY_PATH = pathlib.Path(__file__).parent / 'data_registry.json'
_REPO_ROOT      = pathlib.Path(__file__).parent.parent


# ---------------------------------------------------------------------------
# Registry helpers
# ---------------------------------------------------------------------------

def _load_registry() -> dict:
    with _REGISTRY_PATH.open() as f:
        return json.load(f)


def _save_registry(registry: dict) -> None:
    with _REGISTRY_PATH.open('w') as f:
        json.dump(registry, f, indent=2)
    print(f'  Updated data_registry.json')


def list_datasets() -> list[dict]:
    """Return all registered datasets from data_registry.json."""
    return _load_registry()['datasets']


# ---------------------------------------------------------------------------
# Dataset-specific downloaders
# ---------------------------------------------------------------------------

def _download_global_tectonics(dest: pathlib.Path) -> pathlib.Path:
    """Download global_tectonics from Zenodo and extract to *dest*."""
    zenodo_id  = '6586972'
    api_url    = f'https://zenodo.org/api/records/{zenodo_id}'

    print(f'  Fetching Zenodo record {zenodo_id} …')
    with urllib.request.urlopen(api_url, timeout=30) as r:
        record = json.loads(r.read())

    # Locate the zip (or tar.gz) in the record's file list
    files = record.get('files', [])
    zip_entry = next(
        (f for f in files if f['key'].endswith(('.zip', '.tar.gz'))),
        None,
    )
    if zip_entry is None:
        raise RuntimeError(
            f'No zip/tar.gz found in Zenodo record {zenodo_id}. '
            f'Files available: {[f["key"] for f in files]}'
        )

    download_url = zip_entry['links']['self']
    zip_name     = zip_entry['key']
    zip_path     = dest.parent / zip_name

    print(f'  Downloading {zip_name} ({zip_entry["size"] / 1e6:.1f} MB) …')
    urllib.request.urlretrieve(download_url, zip_path)

    print(f'  Extracting to {dest} …')
    dest.mkdir(parents=True, exist_ok=True)
    if zip_name.endswith('.zip'):
        with zipfile.ZipFile(zip_path, 'r') as z:
            # Strip the top-level directory if the zip has one
            members = z.namelist()
            top_dirs = {m.split('/')[0] for m in members if '/' in m}
            if len(top_dirs) == 1:
                top = top_dirs.pop()
                for member in members:
                    target = dest / pathlib.Path(member).relative_to(top)
                    if member.endswith('/'):
                        target.mkdir(parents=True, exist_ok=True)
                    else:
                        target.parent.mkdir(parents=True, exist_ok=True)
                        with z.open(member) as src, target.open('wb') as dst:
                            shutil.copyfileobj(src, dst)
            else:
                z.extractall(dest)
    else:
        import tarfile
        with tarfile.open(zip_path, 'r:gz') as t:
            t.extractall(dest)

    zip_path.unlink()
    print(f'  global_tectonics installed at {dest}')
    return dest


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

_DOWNLOADERS = {
    'global_tectonics': _download_global_tectonics,
}


def download_dataset(
    name: str,
    dest: str | pathlib.Path | None = None,
) -> pathlib.Path:
    """Download a registered dataset and update the registry with its local path.

    Parameters
    ----------
    name : str
        Dataset name as listed in data_registry.json (e.g.
        ``'global_tectonics'``).
    dest : path-like, optional
        Destination directory.  Defaults to ``<repo_root>/data/<name>/``.

    Returns
    -------
    pathlib.Path
        The local directory / file where the data was installed.

    Raises
    ------
    ValueError
        If *name* is not a registered dataset or has no automated downloader.
    RuntimeError
        If the download or extraction fails.
    """
    registry = _load_registry()
    entry = next((d for d in registry['datasets'] if d['name'] == name), None)
    if entry is None:
        raise ValueError(
            f"Dataset '{name}' not in data_registry.json. "
            f"Known datasets: {[d['name'] for d in registry['datasets']]}"
        )

    if name not in _DOWNLOADERS:
        raise ValueError(
            f"No automated downloader for '{name}'. "
            f"Download it manually — see data_registry.json for the URL/location."
        )

    if dest is None:
        dest = _REPO_ROOT / 'data' / name
    dest = pathlib.Path(dest)

    local_path = _DOWNLOADERS[name](dest)

    # Update the registry so future imports resolve without re-downloading
    entry['location'] = str(local_path)
    _save_registry(registry)

    return local_path
