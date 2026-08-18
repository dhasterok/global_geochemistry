"""
Linear regression tools: ordinary/weighted least squares with confidence
intervals, Deming (errors-in-both-variables) regression, and a convenience
function for plotting a fitted line with its confidence band.
"""

import numpy as np
from scipy import stats


def llsq(y, A, alpha=0.95):
    """General linear least squares fit for a data vector and design matrix.

    Bundles the model parameters, their confidence intervals, RMS misfit,
    and R^2 into a single call, so these diagnostics don't need to be
    re-derived by hand for every fit.

    Parameters
    ----------
    y : array_like, shape (n,)
        Observed data vector.
    A : array_like, shape (n, p)
        Design (operator) matrix, e.g. ``numpy.column_stack([ones, x])``
        for a simple straight-line fit.
    alpha : float, optional
        Confidence level for `ci`, by default 0.95.

    Returns
    -------
    dict
        ``m`` (parameter estimates, shape (p,)), ``ci`` (parameter
        confidence half-widths, shape (p,)), ``rms`` (root-mean-square
        residual), ``rsq`` (coefficient of determination), and
        ``residuals`` (shape (n,)).
    """
    y = np.asarray(y, dtype=float).ravel()
    A = np.asarray(A, dtype=float)
    n, p = A.shape

    C = np.linalg.inv(A.T @ A)
    m = C @ A.T @ y
    y_hat = A @ m
    r = y - y_hat

    df = n - p
    t_c = stats.t.ppf((1 + alpha) / 2, df)
    ci = t_c * np.sqrt(np.diag(C) * (r @ r) / (df - 1))

    rms = np.sqrt(np.sum(r**2) / n)
    rsq = np.sum((y_hat - y.mean())**2) / np.sum((y - y.mean())**2)

    return {'m': m, 'ci': ci, 'rms': float(rms), 'rsq': float(rsq), 'residuals': r}


def linear_fit(x, y, alpha=0.95):
    """Straight-line fit y = m[0] + m[1]*x via :func:`llsq`."""
    x = np.asarray(x, dtype=float).ravel()
    A = np.column_stack([np.ones_like(x), x])
    return llsq(y, A, alpha=alpha)


def wllsq(x, y, sigx, sigy, alpha=0.95, tol=1e-4, max_iter=25):
    """Weighted least squares fit accounting for uncertainty in both x and y.

    Iteratively reweights the fit using the "effective variance" method
    (weights = 1 / (sigma_y^2 + slope^2 * sigma_x^2)), which approximates
    errors-in-both-variables regression without assuming a fixed
    error-variance ratio. For an exact closed-form treatment when the
    ratio is known (or assumed equal), see :func:`deming`.

    Parameters
    ----------
    x, y, sigx, sigy : array_like, shape (n,)
        Data and their 1-sigma uncertainties.
    alpha : float, optional
        Confidence level for `ci`, by default 0.95.
    tol : float, optional
        Relative residual-norm change below which iteration stops.
    max_iter : int, optional
        Maximum number of reweighting iterations.

    Returns
    -------
    dict
        ``m`` (intercept, slope), ``ci`` (confidence half-widths for
        each), ``rsq``, and ``residuals`` (in weighted units).
    """
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    sigx = np.asarray(sigx, dtype=float).ravel()
    sigy = np.asarray(sigy, dtype=float).ravel()

    A = np.column_stack([np.ones_like(x), x])

    w = 1.0 / sigy
    Aw = A * w[:, None]
    yw = y * w
    m = np.linalg.solve(Aw.T @ Aw, Aw.T @ yw)
    r0 = yw - Aw @ m

    C = np.linalg.inv(Aw.T @ Aw)
    for _ in range(max_iter):
        w = 1.0 / np.sqrt(sigy**2 + m[1]**2 * sigx**2)
        Aw = A * w[:, None]
        yw = y * w

        C = np.linalg.inv(Aw.T @ Aw)
        m_new = C @ Aw.T @ yw
        r = yw - Aw @ m_new

        converged = np.linalg.norm(r) / np.linalg.norm(r0) < tol
        m = m_new
        if converged:
            break

    z_c = stats.norm.ppf((1 + alpha) / 2)
    ci = z_c * np.sqrt(np.diag(C))
    rsq = np.sum((Aw @ m - yw.mean())**2) / np.sum((yw - yw.mean())**2)

    return {'m': m, 'ci': ci, 'rsq': float(rsq), 'residuals': r}


def deming(x, y, lam=1.0, ci=False, alpha=0.05):
    """Deming (errors-in-both-variables) linear regression.

    Fits y = b0 + b1*x assuming both x and y contain measurement error,
    related by ``lam = sigma_y^2 / sigma_x^2``.

    *Reference:* Jensen, A.C., 2007, description of the Deming regression
    function for MethComp.

    Parameters
    ----------
    x, y : array_like, shape (n,)
    lam : float, optional
        Ratio of y to x measurement-error variance, by default 1.0.
    ci : bool, optional
        If True, also compute jackknife standard errors and confidence
        intervals for the slope and intercept (an additional n linear fits).
    alpha : float, optional
        Significance level for the confidence interval, used only if `ci`.

    Returns
    -------
    dict
        ``b`` (intercept, slope), ``sigma2_x`` (error variance in x),
        ``cod`` (coefficient of determination), ``x_est``/``y_est``
        (estimated true values), and -- if `ci` -- ``se`` (jackknife
        standard errors) and ``ci`` (confidence half-widths) for (b0, b1).
    """
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    n = len(x)

    cov = np.cov(x, y)
    s_xx, s_xy, s_yy = cov[0, 0], cov[0, 1], cov[1, 1]

    b1 = (s_yy - lam * s_xx + np.sqrt((s_yy - lam * s_xx)**2 + 4 * lam * s_xy**2)) / (2 * s_xy)
    b0 = y.mean() - b1 * x.mean()
    b = np.array([b0, b1])

    x_est = x + b1 / (b1**2 + lam) * (y - b0 - b1 * x)
    y_est = b0 + b1 * x_est

    tmp = np.sum(lam * (x - x_est)**2 + (y - b0 - b1 * x_est)**2)
    sigma2_x = tmp / (2 * lam * (n - 2))
    cod = 1 - tmp / (np.linalg.det(cov) * (n - 1)**2)

    result = {'b': b, 'sigma2_x': float(sigma2_x), 'cod': float(cod), 'x_est': x_est, 'y_est': y_est}

    if ci:
        b_sub = np.empty((n, 2))
        mask = np.ones(n, dtype=bool)
        for i in range(n):
            mask[:] = True
            mask[i] = False
            b_sub[i] = deming(x[mask], y[mask], lam)['b']
        se = b_sub.std(axis=0, ddof=0) * (n - 1) / np.sqrt(n)
        t_c = stats.t.ppf(1 - alpha / 2, n - 2)
        result['se'] = se
        result['ci'] = t_c * se

    return result


def plot_llsq(x, y, m, alpha=0.95, xlim=None, w=None, ax=None, **kwargs):
    """Plots a fitted linear model with a confidence band.

    Given previously fitted coefficients `m` (e.g. from :func:`llsq` or
    :func:`wllsq`), draws the regression line over `xlim` and a shaded
    confidence band derived from the residuals of `y` against `x` -- so
    the fitting method stays decoupled from the plotting.

    Parameters
    ----------
    x, y : array_like, shape (n,)
        Data used to compute the standard error of the fitted line.
    m : array_like, shape (2,) or (2, k)
        Model coefficients (intercept, slope). If 2-D, only the first
        column is used (matches the output of :func:`llsq`/:func:`wllsq`,
        whose second column holds parameter uncertainties, not needed here).
    alpha : float, optional
        Confidence level for the band, by default 0.95.
    xlim : (float, float), optional
        x-range over which to draw the line/band, by default the data range.
    w : array_like, optional
        Per-point weights used only to compute the weighted mean of x for
        the confidence band; by default unweighted.
    ax : matplotlib.axes.Axes, optional
    **kwargs
        Passed to the regression-line plot call.

    Returns
    -------
    matplotlib.axes.Axes
    """
    import matplotlib.pyplot as plt

    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    m = np.asarray(m, dtype=float)
    m0 = m[:, 0] if m.ndim > 1 else m

    w = np.ones_like(x) if w is None else np.asarray(w, dtype=float).ravel()

    if ax is None:
        ax = plt.gca()
    if xlim is None:
        xlim = (x.min(), x.max())

    n_d = len(x)
    n_m = len(m0)
    nu = n_d - n_m

    x_line = np.linspace(xlim[0], xlim[1], 50)
    y_line = m0[0] + m0[1] * x_line

    x_mu = np.sum(w * x) / np.sum(w)
    y_hat = m0[0] + m0[1] * x
    resid_ss = np.sum((y - y_hat)**2)

    se = np.sqrt(resid_ss / nu) * np.sqrt(1 / n_d + (x_line - x_mu)**2 / np.sum((x - x_mu)**2))
    t_c = stats.t.ppf((1 + alpha) / 2, nu)

    ax.fill_between(x_line, y_line - t_c * se, y_line + t_c * se, color='0.7', edgecolor='none')
    ax.plot(x_line, y_line, 'k-', **kwargs)

    return ax
