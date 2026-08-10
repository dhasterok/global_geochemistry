# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# Unlike lame-core/siesta/blueberry, this package's source directory is
# literally named `src` (pyproject.toml maps it to the installed name
# `global_geochemistry` via [tool.setuptools.package-dir]), so a plain
# sys.path insert can't make `import global_geochemistry` resolve here the
# way it does for those. Build this from an environment that has the package
# installed (editable or otherwise) -- LaME's own .venv already does, since
# LaME depends on it: `pip install -e .` from this repo's root if building
# standalone.

project = "global_geochemistry"
copyright = "2025, Derrick Hasterok"
author = "Derrick Hasterok"
release = "0.1.0"

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "numpydoc",
]

autosummary_generate = True
autosummary_imported_members = False
templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**/__pycache__/**", "**/.venv/**", "**/venv/**"]

autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "inherited-members": False,
    "show-inheritance": True,
}
autoclass_content = "both"

# see lame-core's docs/source/conf.py (and LaME's) for why this matters:
# without it, numpydoc's per-class "Methods" table tries to link every
# inherited member to an autosummary stub page that never gets generated,
# producing spurious "stub file not found" warnings.
numpydoc_class_members_toctree = False

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "pandas": ("https://pandas.pydata.org/pandas-docs/stable/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
}

# -- Options for HTML output -------------------------------------------------

html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
html_theme_options = {
    "navbar_start": ["navbar-logo"],
    "navbar_center": ["navbar-nav"],
    "navbar_end": ["navbar-icon-links"],
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/dhasterok/global_geochemistry",
            "icon": "fa-brands fa-github",
            "type": "fontawesome",
        },
    ],
}