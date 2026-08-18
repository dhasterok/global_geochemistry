"""
Compositional data analysis (CoDA).

The log-ratio transforms, perturbation/power operations, and the Aitchison
inner product are provided by the `composition-stats` package (a small,
actively maintained library dedicated to this math) and re-exported here
under this course's established naming convention. `aitchdist` and
`gram_schmidt` are small additions `composition-stats` doesn't provide
directly.

Note: `ilr`/`anti_ilr` here use `composition-stats`' default (Egozcue)
orthonormal basis, which is a valid but not necessarily identical basis to
the one the original MATLAB `ilr.m`/Gram-Schmidt procedure would have
produced -- both preserve all compositional distances/angles, but the
individual ILR coordinate axes may differ in sign or order.

`pfa` and `robust_pca` are working, modern replacements for MATLAB's
`pfa.m`/`robustpca.m`; the original MATLAB sources for both were
incomplete, half-translated ports of R code (`pfa.m` in particular mixes R
and MATLAB syntax throughout and cannot run) rather than functioning
implementations, so these are new implementations built from the stated
intent (ILR + factor analysis, and ILR + a robust covariance estimate)
rather than direct ports.
"""

import numpy as np
from composition_stats import (
    closure,
    clr,
    clr_inv as anti_clr,
    ilr,
    ilr_inv as anti_ilr,
    perturb as aitchperturb,
    power as aitchpower,
    multiplicative_replacement,
)

__all__ = [
    'closure', 'clr', 'anti_clr', 'ilr', 'anti_ilr',
    'aitchperturb', 'aitchpower', 'aitchip', 'aitchdist',
    'multiplicative_replacement', 'gram_schmidt', 'pfa', 'robust_pca',
]


def aitchip(x, y):
    """Aitchison inner product between corresponding rows of x and y.

    ``<x, y>_a = sum(clr(x) * clr(y))``, computed row-wise. Note this
    differs from `composition_stats.inner`, which returns the full
    pairwise ``(n, n)`` dot-product matrix between all rows of `x` and all
    rows of `y` rather than a row-wise ``(n,)`` result -- this wrapper
    matches the row-wise convention of MATLAB's `aitchip.m`
    (``sum(x.*y, 2)``).

    Parameters
    ----------
    x, y : array_like, shape (n, D) or (D,)
        Compositions.

    Returns
    -------
    numpy.ndarray or float
        Aitchison inner product for each row of `x` paired with the
        corresponding row of `y`.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    return np.sum(clr(x) * clr(y), axis=-1)


def aitchdist(x, y=None):
    """Aitchison distance between compositions.

    ``d_a(x, y) = || clr(x) - clr(y) ||``, the Euclidean distance between
    compositions after the centered log-ratio transform.

    Parameters
    ----------
    x : array_like, shape (n, D) or (D,)
        Compositions.
    y : array_like, shape matching x, optional
        Compositions to compare against. Defaults to the uniform (centre)
        composition.

    Returns
    -------
    numpy.ndarray or float
        Aitchison distance for each row of `x`.
    """
    x = np.asarray(x, dtype=float)
    y = np.full_like(x, 1.0 / x.shape[-1]) if y is None else np.asarray(y, dtype=float)
    return np.linalg.norm(clr(x) - clr(y), axis=-1)


def gram_schmidt(A):
    """Gram-Schmidt orthonormalization of the columns of A.

    Not needed for :func:`ilr`/:func:`anti_ilr` (which build their own
    default orthonormal basis internally), but included as a standalone
    teaching utility for constructing custom log-ratio bases by hand.

    Parameters
    ----------
    A : array_like, shape (m, n)

    Returns
    -------
    Q : numpy.ndarray, shape (m, n)
        Orthonormal basis vectors (columns).
    R : numpy.ndarray, shape (n, n)
        Upper-triangular coefficients such that ``A = Q @ R``.
    """
    A = np.asarray(A, dtype=float)
    m, n = A.shape
    Q = np.zeros((m, n))
    R = np.zeros((n, n))
    for j in range(n):
        v = A[:, j].copy()
        for i in range(j):
            R[i, j] = Q[:, i] @ A[:, j]
            v -= R[i, j] * Q[:, i]
        R[j, j] = np.linalg.norm(v)
        Q[:, j] = v / R[j, j]
    return Q, R


def _prepare(data, columns):
    import pandas as pd
    if isinstance(data, pd.DataFrame):
        columns = list(data.columns) if columns is None else columns
        x = data[columns].to_numpy(dtype=float)
    else:
        x = np.asarray(data, dtype=float)
    return closure(x)


def pfa(data, columns=None, n_factors=2):
    """Principal factor analysis of compositional data.

    Transforms the data via the isometric log-ratio (ILR) transform, then
    fits a maximum-likelihood factor model
    (`sklearn.decomposition.FactorAnalysis`) -- a working, modern
    equivalent of MATLAB's `pfa.m`, whose source was an incomplete,
    non-functional port of R's `factanal` (see module docstring).

    Parameters
    ----------
    data : pandas.DataFrame or array_like
        Compositional data (need not already be closed to a constant sum).
    columns : list of str, optional
        Column names to use, if `data` is a DataFrame.
    n_factors : int, optional
        Number of factors to extract, by default 2.

    Returns
    -------
    dict
        ``loadings`` (shape (n_factors, n_features - 1), in ILR space),
        ``scores`` (shape (n_samples, n_factors)), ``noise_variance``,
        and ``z`` (the ILR-transformed data used).
    """
    from sklearn.decomposition import FactorAnalysis

    z = ilr(_prepare(data, columns))
    model = FactorAnalysis(n_components=n_factors)
    scores = model.fit_transform(z)

    return {
        'loadings': model.components_,
        'scores': scores,
        'noise_variance': model.noise_variance_,
        'z': z,
    }


def robust_pca(data, columns=None, n_components=None):
    """Robust PCA on compositional data via the ILR transform.

    Transforms the data via the ILR transform, estimates a robust
    covariance and location with the Minimum Covariance Determinant
    estimator (`sklearn.covariance.MinCovDet`), then extracts principal
    components by eigendecomposition of that robust covariance -- so a
    handful of outlying compositions cannot dominate the components. This
    is a working, modern equivalent of MATLAB's `robustpca.m` (ILR + MCD
    covariance + `pfa`), whose source called into the same broken `pfa.m`
    described in the module docstring.

    Parameters
    ----------
    data : pandas.DataFrame or array_like
        Compositional data (need not already be closed to a constant sum).
    columns : list of str, optional
        Column names to use, if `data` is a DataFrame.
    n_components : int, optional
        Number of components to keep. Defaults to all (n_features - 1
        after the ILR transform).

    Returns
    -------
    dict
        ``scores``, ``components`` (shape (n_components, n_features - 1),
        matching `sklearn.decomposition.PCA`'s convention),
        ``explained_variance``, ``explained_variance_ratio``,
        ``location`` (robust center in ILR space), and ``z`` (the
        ILR-transformed data used).
    """
    from sklearn.covariance import MinCovDet

    z = ilr(_prepare(data, columns))

    mcd = MinCovDet().fit(z)
    eigvals, eigvecs = np.linalg.eigh(mcd.covariance_)
    order = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]

    n_components = n_components or z.shape[1]
    components = eigvecs[:, :n_components]
    explained_variance = eigvals[:n_components]

    return {
        'scores': (z - mcd.location_) @ components,
        'components': components.T,
        'explained_variance': explained_variance,
        'explained_variance_ratio': explained_variance / eigvals.sum(),
        'location': mcd.location_,
        'z': z,
    }
