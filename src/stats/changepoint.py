import numpy as np
from scipy import stats


def pettitt(x):
    """Pettitt's test for a single change point in a univariate series.

    A nonparametric, rank-based test of the null hypothesis that a series
    has no change point against the alternative that it does. Reports only
    the single most likely change point and gives no indication of how
    probable a change is at other locations -- for multiple change points or
    a full probability profile, use :func:`ruptures_changepoints` or
    :func:`bocpd_changepoints` instead.

    *Reference:* Pohlert, T., 2016, Non-Parametric Trend Tests and
    Change-Point Detection.

    Parameters
    ----------
    x : array_like
        Univariate data series.

    Returns
    -------
    dict
        ``location`` (index of the change point, or ``None`` if the test
        statistic is tied across more than one location), ``statistic``
        (Pettitt's K statistic), and ``p_value``.
    """
    x = np.asarray(x, dtype=float).ravel()
    m = len(x)

    V = np.sign(x[:, None] - x[None, :]).sum(axis=1)
    U = np.cumsum(V)[:-1]

    K = np.max(np.abs(U))
    loc = np.flatnonzero(np.abs(U) == K)
    p_value = 2 * np.exp(-6 * K**2 / (m**3 + m**2))

    if len(loc) > 1:
        return {'location': None, 'statistic': np.nan, 'p_value': np.nan}
    return {'location': int(loc[0]), 'statistic': float(K), 'p_value': float(p_value)}


def ruptures_changepoints(x, method='pelt', model='rbf', pen=10, n_bkps=None, min_size=2, jump=1):
    """Detects multiple change points using the `ruptures` package.

    Unlike :func:`pettitt`, this can return more than one change point, at
    the cost of a deterministic (non-probabilistic) location estimate.

    Parameters
    ----------
    x : array_like
        Univariate or multivariate data series (samples along axis 0).
    method : {'pelt', 'binseg', 'bottomup', 'window', 'dynp'}, optional
        `ruptures` search algorithm, by default 'pelt'.
    model : str, optional
        Cost function name understood by `ruptures` (e.g. 'rbf', 'l1',
        'l2'), by default 'rbf'.
    pen : float, optional
        Penalty used by penalty-based methods ('pelt', 'binseg', 'bottomup',
        'window') to control the number of change points found. Ignored
        when `n_bkps` is given.
    n_bkps : int, optional
        Exact number of change points to find. Required for 'dynp'; for
        other methods this overrides `pen`.
    min_size, jump : int, optional
        Passed through to the `ruptures` algorithm constructor.

    Returns
    -------
    list of int
        Indices of the detected change points.
    """
    import ruptures as rpt

    x = np.asarray(x, dtype=float)
    signal = x.reshape(-1, 1) if x.ndim == 1 else x

    algorithms = {
        'pelt': rpt.Pelt,
        'binseg': rpt.Binseg,
        'bottomup': rpt.BottomUp,
        'window': rpt.Window,
        'dynp': rpt.Dynp,
    }
    if method not in algorithms:
        raise ValueError(f"Unknown method: {method!r}. Choose from {list(algorithms)}.")

    algo = algorithms[method](model=model, min_size=min_size, jump=jump).fit(signal)

    if method == 'dynp' or n_bkps is not None:
        if n_bkps is None:
            raise ValueError("n_bkps is required for method='dynp'.")
        breakpoints = algo.predict(n_bkps=n_bkps)
    else:
        breakpoints = algo.predict(pen=pen)

    # ruptures appends len(signal) as a trailing sentinel marking the end of
    # the series; it is not a change point
    return [b for b in breakpoints if b < len(x)]


def bocpd(x, hazard=1 / 250, mu0=None, kappa0=1.0, alpha0=1.0, beta0=1.0):
    """Bayesian online change point detection for a univariate series.

    Implements the algorithm of Adams & MacKay (2007), assuming a Gaussian
    generative model with unknown mean and variance (Normal-Gamma conjugate
    prior) and a constant hazard rate (equivalent to a geometric prior on
    run length). Unlike :func:`pettitt`, this yields a probability of a
    change point at *every* location, from which multiple change points can
    be read off -- see :func:`bocpd_changepoints`.

    *Reference:* Adams, R.P., MacKay, D.J.C., 2007, Bayesian Online
    Changepoint Detection, arXiv:0710.3742.

    Parameters
    ----------
    x : array_like
        Univariate observation series.
    hazard : float, optional
        Constant hazard rate, i.e. the prior probability of a change point
        at any given time step. Equivalent to 1 / (expected run length).
    mu0, kappa0, alpha0, beta0 : float, optional
        Normal-Gamma prior hyperparameters for the unknown mean/precision
        of each segment. `mu0` defaults to `x[0]`.

    Returns
    -------
    cp_prob : numpy.ndarray, shape (len(x),)
        ``P(r_t = 0 | x_{1:t})`` at each t, as defined in Adams & MacKay's
        Algorithm 1. Note that for a *constant* hazard rate this quantity is
        a mathematical identity equal to `hazard` at every t, regardless of
        the data -- it is returned for completeness/comparison with the
        paper, but is not by itself useful for locating change points. Use
        `run_length_posterior` (or :func:`bocpd_changepoints`) instead.
    run_length_posterior : numpy.ndarray, shape (len(x), len(x))
        Full run-length posterior distribution; ``run_length_posterior[t, r]``
        is the probability that the run length at time `t` is `r`. This is
        the informative output -- e.g. its row-wise argmax is the maximum a
        posteriori (MAP) run length trajectory, which resets sharply at
        genuine change points.
    """
    x = np.asarray(x, dtype=float).ravel()
    T = len(x)
    if mu0 is None:
        mu0 = x[0]

    mu = np.array([mu0])
    kappa = np.array([kappa0])
    alpha = np.array([alpha0])
    beta = np.array([beta0])

    R = np.zeros((T + 1, T + 1))
    R[0, 0] = 1.0
    cp_prob = np.zeros(T)

    for t in range(T):
        df = 2 * alpha
        scale = np.sqrt(beta * (kappa + 1) / (alpha * kappa))
        pred = stats.t.pdf(x[t], df, loc=mu, scale=scale)

        growth = R[t, :t + 1] * pred * (1 - hazard)
        cp = np.sum(R[t, :t + 1] * pred * hazard)

        R[t + 1, 1:t + 2] = growth
        R[t + 1, 0] = cp
        R[t + 1, :t + 2] /= R[t + 1, :t + 2].sum()

        cp_prob[t] = R[t + 1, 0]

        new_kappa = np.concatenate(([kappa0], kappa + 1))
        new_mu = np.concatenate(([mu0], (kappa * mu + x[t]) / (kappa + 1)))
        new_alpha = np.concatenate(([alpha0], alpha + 0.5))
        new_beta = np.concatenate((
            [beta0],
            beta + kappa * (x[t] - mu) ** 2 / (2 * (kappa + 1)),
        ))
        mu, kappa, alpha, beta = new_mu, new_kappa, new_alpha, new_beta

    return cp_prob, R[1:, :T]


def bocpd_changepoints(x, hazard=1 / 250, window=3, min_probability=0.5, **kwargs):
    """Runs :func:`bocpd` and extracts discrete change points with probabilities.

    The marginal ``P(r_t = 0 | x_{1:t})`` returned by :func:`bocpd` is fixed
    at exactly `hazard` for every t under a constant hazard rate -- a
    mathematical property of the model, not a data-driven signal -- so it
    cannot be used to locate change points on its own. Instead, change
    points are identified from where the *most probable* run length
    (``argmax_r run_length_posterior[t, :]``) fails to simply grow by one
    step, i.e. the model's belief genuinely restarts. Each detected reset is
    assigned a probability equal to the posterior mass on a short run
    length, ``P(r_t <= window | x_{1:t})``, which does vary meaningfully
    with the data and is high exactly at genuine breaks.

    Parameters
    ----------
    x : array_like
        Univariate observation series.
    hazard : float, optional
        See :func:`bocpd`.
    window : int, optional
        Number of short run lengths (``0..window``) aggregated into the
        reported probability.
    min_probability : float, optional
        Minimum ``P(r_t <= window)`` for a detected reset to be reported.
    **kwargs
        Additional keyword arguments passed to :func:`bocpd` (`mu0`,
        `kappa0`, `alpha0`, `beta0`).

    Returns
    -------
    list of dict
        Each entry has ``location`` (index into `x`) and ``probability``.
    """
    _, R = bocpd(x, hazard=hazard, **kwargs)
    map_run_length = np.argmax(R, axis=1)
    p_short_run = R[:, :window + 1].sum(axis=1)

    resets = np.flatnonzero(np.diff(map_run_length, prepend=map_run_length[0] - 1) <= 0)
    return [
        {'location': int(t), 'probability': float(p_short_run[t])}
        for t in resets if p_short_run[t] >= min_probability
    ]
