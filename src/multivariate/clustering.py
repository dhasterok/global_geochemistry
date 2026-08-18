"""
Hierarchical clustering for geochemical (or other numeric) data.

Generalizes the pattern in the MATLAB `dendrite` function -- cluster
samples, then examine how cluster membership relates to an external
group label such as rock type -- without hard-coding a specific dataset's
columns.
"""

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.spatial.distance import pdist


def hierarchical_clustering(data, method='average', metric='euclidean', n_clusters=None):
    """Agglomerative hierarchical clustering via `scipy.cluster.hierarchy`.

    Parameters
    ----------
    data : array_like, shape (n_samples, n_features)
    method : str, optional
        Linkage method (e.g. 'average', 'single', 'complete', 'ward'),
        by default 'average'.
    metric : str, optional
        Distance metric passed to `scipy.spatial.distance.pdist`, by
        default 'euclidean'.
    n_clusters : int, optional
        If given, also returns flat cluster assignments cut at this many
        clusters (via `scipy.cluster.hierarchy.fcluster`).

    Returns
    -------
    dict
        ``linkage`` (the linkage matrix Z) and, if `n_clusters` is given,
        ``clusters`` (1-based cluster id per sample, shape (n_samples,)).
    """
    x = np.asarray(data, dtype=float)
    y = pdist(x, metric=metric)
    z = linkage(y, method=method)

    result = {'linkage': z}
    if n_clusters is not None:
        result['clusters'] = fcluster(z, n_clusters, criterion='maxclust')
    return result


def plot_dendrogram(linkage_matrix, ax=None, **kwargs):
    """Plots a dendrogram from a linkage matrix.

    Thin convenience wrapper over `scipy.cluster.hierarchy.dendrogram` so
    students don't need to import it separately after
    :func:`hierarchical_clustering`.

    Parameters
    ----------
    linkage_matrix : array_like or dict
        A linkage matrix Z, or the dict returned by
        :func:`hierarchical_clustering` (uses its ``linkage`` entry).
    ax : matplotlib.axes.Axes, optional
    **kwargs
        Passed to `scipy.cluster.hierarchy.dendrogram`.

    Returns
    -------
    matplotlib.axes.Axes
    dict
        The dendrogram layout data returned by `scipy`'s `dendrogram`.
    """
    import matplotlib.pyplot as plt

    z = linkage_matrix['linkage'] if isinstance(linkage_matrix, dict) else linkage_matrix

    if ax is None:
        ax = plt.gca()
    ddata = dendrogram(z, ax=ax, **kwargs)
    return ax, ddata


def cluster_group_composition(clusters, groups, ax=None, normalize=True, cmap='viridis'):
    """Plots each group's distribution across clusters as a heatmap.

    E.g. with `groups` as rock type and `clusters` from
    :func:`hierarchical_clustering`, shows what fraction of each rock
    type's samples fall into each cluster -- the generalized form of the
    MATLAB `dendrite` function's rock-type-vs-cluster panel.

    Parameters
    ----------
    clusters : array_like, shape (n,)
        Cluster assignment per sample.
    groups : array_like, shape (n,)
        Group label per sample (e.g. rock type).
    ax : matplotlib.axes.Axes, optional
    normalize : bool, optional
        If True (default), normalize each group's row to sum to 1.
    cmap : str, optional

    Returns
    -------
    matplotlib.axes.Axes
    pandas.DataFrame
        Group-by-cluster count (or fraction) table.
    """
    import matplotlib.pyplot as plt

    clusters = np.asarray(clusters)
    groups = np.asarray(groups)

    cluster_ids = np.unique(clusters)
    group_ids = np.unique(groups)

    table = pd.DataFrame(0.0, index=group_ids, columns=cluster_ids)
    for g in group_ids:
        for c in cluster_ids:
            table.loc[g, c] = np.sum((groups == g) & (clusters == c))

    if normalize:
        table = table.div(table.sum(axis=1), axis=0)

    if ax is None:
        _, ax = plt.subplots()

    im = ax.imshow(table.to_numpy(), aspect='auto', cmap=cmap)
    ax.set_xticks(range(len(cluster_ids)))
    ax.set_xticklabels(cluster_ids)
    ax.set_yticks(range(len(group_ids)))
    ax.set_yticklabels(group_ids)
    ax.set_xlabel('Cluster')
    plt.colorbar(im, ax=ax)

    return ax, table
