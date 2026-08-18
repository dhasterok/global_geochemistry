"""
Principal component analysis for geochemical (or other numeric) data.

`pca_variance` and `pca_vectors` are adapted from LaserMapExplorer's
`plot_pca_variance` / `plot_pca_components` (src/plotting/LamePlot.py),
generalized to take a plain fitted `sklearn.decomposition.PCA` object
instead of the app's style/canvas/data-handler state, so the same
implementation can be shared between the two projects.
"""

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def pca(data, columns=None, n_components=None, scale=True, **pca_kwargs):
    """Runs PCA on a table of numeric data.

    Standardizes each column to zero mean/unit variance before fitting
    (the common convention when columns are in different units or have
    very different scales), unless `scale=False`.

    Parameters
    ----------
    data : pandas.DataFrame or array_like
        Numeric data, samples along rows.
    columns : list of str, optional
        Column names to use, if `data` is a DataFrame. Defaults to all
        columns.
    n_components : int, optional
        Number of components to keep. Defaults to
        ``min(n_samples, n_features)``.
    scale : bool, optional
        If True (default), standardize columns via
        `sklearn.preprocessing.StandardScaler` before fitting.
    **pca_kwargs
        Additional keyword arguments passed to `sklearn.decomposition.PCA`,
        e.g. ``svd_solver='randomized'`` for large datasets (equivalent to
        the MATLAB `rPCA` randomized approximation).

    Returns
    -------
    dict
        ``pca`` (the fitted `sklearn.decomposition.PCA` object), ``scores``
        (ndarray, shape (n_samples, n_components)), and ``labels`` (the
        column names/labels used).
    """
    if isinstance(data, pd.DataFrame):
        columns = list(data.columns) if columns is None else columns
        x = data[columns].to_numpy(dtype=float)
    else:
        x = np.asarray(data, dtype=float)
        columns = columns or [str(i) for i in range(x.shape[1])]

    if scale:
        x = StandardScaler().fit_transform(x)

    if n_components is None:
        n_components = min(x.shape)

    model = PCA(n_components=n_components, **pca_kwargs)
    scores = model.fit_transform(x)

    return {'pca': model, 'scores': scores, 'labels': columns}


def _get_model(pca_result):
    return pca_result['pca'] if isinstance(pca_result, dict) else pca_result


def pca_variance(pca_result, ax=None, **kwargs):
    """Plots individual and cumulative explained variance vs. component number.

    Adapted from LaserMapExplorer's `plot_pca_variance` (a scree plot).

    Parameters
    ----------
    pca_result : sklearn.decomposition.PCA or dict
        A fitted PCA object, or the dict returned by :func:`pca`.
    ax : matplotlib.axes.Axes, optional
    **kwargs
        Passed to both line-plot calls (e.g. `linewidth`).

    Returns
    -------
    matplotlib.axes.Axes
    """
    import matplotlib.pyplot as plt

    model = _get_model(pca_result)

    if ax is None:
        ax = plt.gca()

    variances = model.explained_variance_ratio_
    n_components = np.arange(1, len(variances) + 1)
    cumulative = variances.cumsum()

    ax.plot(n_components, variances, linestyle='-', marker='o', label='Individual', **kwargs)
    ax.plot(n_components, cumulative, linestyle='-', marker='s', label='Cumulative', **kwargs)

    ax.set_ylim(0, 1.0)
    ax.set_xlabel('Principal Component')
    ax.set_ylabel('Explained Variance Ratio')
    ax.legend()
    ax.set_box_aspect(1)

    return ax


def pca_vectors(pca_result, labels=None, pc_x=0, pc_y=1, scale=None, ax=None, color='0.3', fontsize=8):
    """Overlays PCA loading vectors (as labeled arrows) on the current axes.

    Adapted from LaserMapExplorer's `plot_pca_components`, generalized to
    take a fitted `sklearn.decomposition.PCA` object and plain variable
    labels rather than the app's data/state objects. Typically drawn on
    top of a score scatter from :func:`pca_score_scatter` to form a biplot.

    Parameters
    ----------
    pca_result : sklearn.decomposition.PCA or dict
        A fitted PCA object, or the dict returned by :func:`pca`.
    labels : list of str, optional
        Variable names, one per column of the original data. Defaults to
        the ``labels`` entry of `pca_result` if it's a dict, else integer
        indices.
    pc_x, pc_y : int, optional
        Zero-based component indices for the two axes, by default (0, 1),
        i.e. PC1 vs PC2.
    scale : float, optional
        Arrow-length multiplier. Defaults to a length that keeps the
        longest arrow within the current axes limits.
    ax : matplotlib.axes.Axes, optional
    color : color spec, optional
    fontsize : int, optional

    Returns
    -------
    matplotlib.axes.Axes
    """
    import matplotlib.pyplot as plt

    model = _get_model(pca_result)
    if labels is None and isinstance(pca_result, dict):
        labels = pca_result.get('labels')

    components = model.components_
    n_vars = components.shape[1]
    if labels is None:
        labels = [str(i) for i in range(n_vars)]

    x = components[pc_x, :]
    y = components[pc_y, :]

    if ax is None:
        ax = plt.gca()

    if scale is None:
        xlim, ylim = ax.get_xlim(), ax.get_ylim()
        max_extent = max(abs(v) for v in (*xlim, *ylim))
        max_len = np.max(np.sqrt(x**2 + y**2)) or 1.0
        scale = 0.8 * max_extent / max_len

    ax.quiver(np.zeros(n_vars), np.zeros(n_vars), scale * x, scale * y, color=color,
              angles='xy', scale_units='xy', scale=1, headlength=3, headaxislength=3)

    for i, label in enumerate(labels):
        ha = 'left' if x[i] >= 0 else 'right'
        va = 'bottom' if y[i] >= 0 else 'top'
        ax.text(scale * x[i], scale * y[i], label, fontsize=fontsize, ha=ha, va=va, color=color)

    return ax


def pca_score_scatter(pca_result, pc_x=0, pc_y=1, groups=None, average=False, ax=None, **kwargs):
    """Scatter plot of PCA scores for two components, optionally grouped.

    Parameters
    ----------
    pca_result : dict or array_like
        The dict returned by :func:`pca` (uses its ``scores``), or a raw
        (n_samples, n_components) score array.
    pc_x, pc_y : int, optional
        Zero-based component indices for the two axes, by default (0, 1).
    groups : array_like, optional
        Per-sample group label used to color points (or compute averages).
    average : bool, optional
        If True and `groups` is given, plot one averaged point per group
        instead of every sample.
    ax : matplotlib.axes.Axes, optional
    **kwargs
        Passed to `Axes.scatter`.

    Returns
    -------
    matplotlib.axes.Axes
    """
    import matplotlib.pyplot as plt

    scores = pca_result['scores'] if isinstance(pca_result, dict) else np.asarray(pca_result)
    x = scores[:, pc_x]
    y = scores[:, pc_y]

    if ax is None:
        ax = plt.gca()

    kwargs.setdefault('alpha', 0.3)
    kwargs.setdefault('s', 12)

    if groups is None:
        ax.scatter(x, y, **kwargs)
    else:
        groups = np.asarray(groups)
        for g in np.unique(groups):
            mask = groups == g
            if average:
                ax.scatter(x[mask].mean(), y[mask].mean(), label=str(g), **kwargs)
            else:
                ax.scatter(x[mask], y[mask], label=str(g), **kwargs)
        ax.legend()

    ax.set_xlabel(f'PC{pc_x + 1}')
    ax.set_ylabel(f'PC{pc_y + 1}')
    ax.set_box_aspect(1)

    return ax
