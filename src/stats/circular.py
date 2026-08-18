import numpy as np


def _u2_statistic(a, b):
    """Watson's U^2 two-sample test statistic for a pair of angular samples."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    ni, nj = len(a), len(b)
    n = ni + nj

    values = np.concatenate([a, b])
    labels = np.concatenate([np.ones(ni), np.zeros(nj)])

    order = np.argsort(values, kind='mergesort')
    values = values[order]
    labels = labels[order]

    ci = np.cumsum(labels)
    cj = np.cumsum(1 - labels)

    # collapse tied angles to the last occurrence of each unique value
    _, last_idx = np.unique(values[::-1], return_index=True)
    last_idx = np.sort(len(values) - 1 - last_idx)

    fi = ci[last_idx] / ni
    fj = cj[last_idx] / nj
    d = fi - fj

    return ni * nj / n**2 * (np.sum(d**2) - np.sum(d)**2 / n)


def watsons_u2(a, b, alpha=None, n_permutations=1000, rng=None):
    """Watson's U^2 two-sample test for circular (angular) data.

    Tests whether two samples of angular data (in the same units, e.g.
    degrees or radians) are drawn from the same circular distribution. The
    null distribution is built by permutation rather than by interpolating
    a fixed critical-value table, so it adapts to the actual sample sizes.

    Parameters
    ----------
    a, b : array_like
        Angular samples for the two groups, in the same units.
    alpha : float, optional
        If given, return the critical U^2 value at significance level
        `alpha` instead of a p-value.
    n_permutations : int, optional
        Number of permutations used to build the null distribution.
    rng : numpy.random.Generator, optional
        Random generator to use, by default a fresh `numpy.random.default_rng()`.

    Returns
    -------
    float
        Observed U^2 statistic.
    float
        p-value (if `alpha` is None) or the critical U^2 value at `alpha`.
    """
    if rng is None:
        rng = np.random.default_rng()

    u2_obs = _u2_statistic(a, b)

    combined = np.concatenate([np.asarray(a, dtype=float), np.asarray(b, dtype=float)])
    n_a = len(a)
    null = np.empty(n_permutations)
    for i in range(n_permutations):
        perm = rng.permutation(combined)
        null[i] = _u2_statistic(perm[:n_a], perm[n_a:])

    if alpha is not None:
        return u2_obs, float(np.quantile(null, 1 - alpha))
    return u2_obs, float(np.mean(null >= u2_obs))
