import numpy as np


def mahalanobis(x, mu=None, cov=None, robust=False, recenter=None):
    """Computes the squared Mahalanobis distance of each row of x to a centroid.

    Parameters
    ----------
    x : array_like, shape (n_samples, n_features)
        Data matrix.
    mu : array_like, optional
        Centroid to measure distance to. Defaults to the (robust or
        classical) sample mean of `x`, after any `recenter` transform.
    cov : array_like, optional
        Covariance matrix. Defaults to the (robust or classical) sample
        covariance of `x`, after any `recenter` transform.
    robust : bool, optional
        If True and `mu`/`cov` are not supplied, estimate them with the
        Minimum Covariance Determinant estimator
        (`sklearn.covariance.MinCovDet`) instead of the classical sample
        mean/covariance, so a handful of outliers cannot dominate the
        distance metric used to find them.
    recenter : {'clr', 'ilr'}, optional
        If given, apply the named compositional log-ratio transform (see
        `geochem.coda`) to `x` before computing distances.

    Returns
    -------
    d_squared : numpy.ndarray, shape (n_samples,)
        Squared Mahalanobis distance for each row of `x`.
    mu : numpy.ndarray
        Centroid used.
    cov : numpy.ndarray
        Covariance matrix used.
    """
    x = np.asarray(x, dtype=float)

    if recenter is not None:
        if recenter == 'clr':
            from geochem.coda import clr
            x = clr(x)
        elif recenter == 'ilr':
            from geochem.coda import ilr
            x = ilr(x)
        else:
            raise ValueError(f"Unknown recenter transform: {recenter!r}")

    if mu is None or cov is None:
        if robust:
            from sklearn.covariance import MinCovDet
            estimator = MinCovDet().fit(x)
            if mu is None:
                mu = estimator.location_
            if cov is None:
                cov = estimator.covariance_
        else:
            if mu is None:
                mu = x.mean(axis=0)
            if cov is None:
                cov = np.cov(x, rowvar=False)

    mu = np.asarray(mu, dtype=float)
    cov = np.atleast_2d(np.asarray(cov, dtype=float))

    delta = x - mu
    inv_cov = np.linalg.inv(cov)
    d_squared = np.einsum('ij,jk,ik->i', delta, inv_cov, delta)

    return d_squared, mu, cov
