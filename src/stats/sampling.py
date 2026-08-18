import numpy as np


def rejsample(x, p, n, kind=None, rng=None):
    """Draws random samples from an arbitrary PDF using rejection sampling.

    Parameters
    ----------
    x : array_like
        Support points, strictly increasing. Interpretation depends on ``kind``:
        for ``'smooth'``, ``x`` gives the locations at which ``p`` is evaluated;
        for ``'hist'``, ``x`` gives the bin edges of a histogram with ``len(p)`` bins.
    p : array_like
        PDF values at ``x`` (``'smooth'``) or histogram bin densities (``'hist'``).
    n : int
        Number of samples to draw.
    kind : {'smooth', 'hist'}, optional
        Distribution representation. Inferred from the relative lengths of ``x``
        and ``p`` when not given: ``len(x) == len(p)`` implies ``'smooth'``,
        ``len(x) == len(p) + 1`` implies ``'hist'``.
    rng : numpy.random.Generator, optional
        Random generator to use, by default a fresh ``numpy.random.default_rng()``.

    Returns
    -------
    numpy.ndarray
        Array of ``n`` samples drawn from the distribution defined by ``(x, p)``.
    """
    x = np.asarray(x, dtype=float)
    p = np.asarray(p, dtype=float)

    if kind is None:
        if len(x) == len(p):
            kind = 'smooth'
        elif len(x) == len(p) + 1:
            kind = 'hist'
        else:
            raise ValueError("x and p must have equal length ('smooth') or x one longer than p ('hist').")

    if np.any(np.diff(x) <= 0):
        raise ValueError('x must be strictly increasing.')

    if rng is None:
        rng = np.random.default_rng()

    if kind == 'smooth':
        return _smooth_sample(x, p, n, rng)
    elif kind == 'hist':
        return _hist_sample(x, p, n, rng)
    else:
        raise ValueError(f"Unknown distribution type: {kind!r}")


def _smooth_sample(x, p, n, rng):
    # resample onto a uniform grid at the finest spacing present, since the
    # envelope/interpolation below assumes evenly spaced support points
    dx = np.diff(x).min()
    xu = np.arange(x[0], x[-1] + dx / 2, dx)
    pu = np.interp(xu, x, p)

    lo, hi = xu[0], xu[-1]
    M = pu.max()

    samples = np.empty(n)
    filled = 0
    while filled < n:
        batch = n - filled
        xr = rng.uniform(lo, hi, size=batch)
        yr = rng.uniform(0, M, size=batch)
        f = np.interp(xr, xu, pu)
        accepted = xr[yr < f]
        k = len(accepted)
        samples[filled:filled + k] = accepted
        filled += k
    return samples


def _hist_sample(x, p, n, rng):
    lo, hi = x[0], x[-1]
    M = p.max()

    samples = np.empty(n)
    filled = 0
    while filled < n:
        batch = n - filled
        xr = rng.uniform(lo, hi, size=batch)
        yr = rng.uniform(0, M, size=batch)
        bin_idx = np.clip(np.searchsorted(x, xr, side='right') - 1, 0, len(p) - 1)
        accepted = xr[yr < p[bin_idx]]
        k = len(accepted)
        samples[filled:filled + k] = accepted
        filled += k
    return samples
