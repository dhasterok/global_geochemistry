import numpy as np
from scipy import stats


def robust_corrcoef(x, y):
    """Computes four correlation coefficients for a pair of variables.

    Bundles Blomqvist's quadrant correlation coefficient (a simple,
    outlier-resistant measure based only on the sign of each point's
    deviation from the median) alongside the standard Spearman, Pearson,
    and Kendall coefficients, so the four can be compared directly.

    Parameters
    ----------
    x, y : array_like
        Equal-length data series.

    Returns
    -------
    dict
        Keys ``quadrant``, ``spearman``, ``pearson``, ``kendall``.
    """
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    n = len(x)

    xmed = np.median(x)
    ymed = np.median(y)
    rq = np.sum(np.sign(x - xmed) * np.sign(y - ymed)) / n

    rs, _ = stats.spearmanr(x, y)
    rp, _ = stats.pearsonr(x, y)
    rk, _ = stats.kendalltau(x, y)

    return {'quadrant': float(rq), 'spearman': float(rs), 'pearson': float(rp), 'kendall': float(rk)}
