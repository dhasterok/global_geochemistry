"""
Geochemical forward models for nonlinear inversion.

Each model function has the signature ``(x, *params) -> y`` so it can
be passed directly to :func:`scipy.optimize.curve_fit`.  A companion
``_jac`` variant provides the Jacobian for use with
:func:`scipy.optimize.least_squares` (``jac=`` keyword).

Available models
----------------
``erf_forward`` / ``erf_jac``
    Error-function model:  y = A · erf(x / w)

``lognormal_forward`` / ``lognormal_jac``
    Log-normal probability density:
    y = exp(-(ln x - μ)² / (2σ²)) / (x √(2π) σ)

Convenience fitting wrappers
-----------------------------
``fit_erf(x, y, p0=None)``
    Fit the error-function model via :func:`scipy.optimize.curve_fit`.

``fit_lognormal(x, y, p0=None)``
    Fit the log-normal model via :func:`scipy.optimize.curve_fit`.
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erf
from scipy.stats import norm, maxwell


# ---------------------------------------------------------------------------
# Error-function model
# ---------------------------------------------------------------------------

def erf_forward(x, A, w):
    """Error-function forward model: y = A · erf(x / w)."""
    return A * erf(x / w)


def erf_jac(x, A, w):
    """Jacobian of erf_forward w.r.t. (A, w).

    Returns an (N, 2) array: [∂y/∂A, ∂y/∂w].
    """
    x = np.asarray(x, dtype=float)
    z = x / w
    dydA = erf(z)
    dydw = -A * (2 / np.sqrt(np.pi)) * np.exp(-z**2) * x / w**2
    return np.column_stack([dydA, dydw])


def fit_erf(x, y, p0=None, **kwargs):
    """Fit error-function model A · erf(x/w) to (x, y) data.

    Parameters
    ----------
    x, y : array-like
    p0 : (A0, w0), optional
        Initial parameter guess.  Defaults to ``(max(y), std(x))``.
    **kwargs
        Forwarded to :func:`scipy.optimize.curve_fit`.

    Returns
    -------
    popt : ndarray, shape (2,)
        Optimal (A, w).
    pcov : ndarray, shape (2, 2)
        Covariance of popt.
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    if p0 is None:
        p0 = [np.nanmax(np.abs(y)), np.nanstd(x) or 1.0]
    return curve_fit(erf_forward, x, y, p0=p0, **kwargs)


# ---------------------------------------------------------------------------
# Log-normal probability density model
# ---------------------------------------------------------------------------

def lognormal_forward(x, mu, sigma):
    """Log-normal PDF forward model.

    y = exp(-(ln(x) - μ)² / (2σ²)) / (x √(2π) σ)
    """
    x = np.asarray(x, dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        ln_x = np.log(x)
        a = (ln_x - mu) / sigma
        return np.exp(-a**2 / 2) / (x * np.sqrt(2 * np.pi) * sigma)


def lognormal_jac(x, mu, sigma):
    """Jacobian of lognormal_forward w.r.t. (μ, σ).

    Returns an (N, 2) array: [∂y/∂μ, ∂y/∂σ].
    """
    x = np.asarray(x, dtype=float)
    y = lognormal_forward(x, mu, sigma)
    a = (np.log(x) - mu) / sigma
    dydmu   = y * a / sigma
    dydsigma = y * (a**2 - 1) / sigma
    return np.column_stack([dydmu, dydsigma])


def fit_lognormal(x, y, p0=None, **kwargs):
    """Fit log-normal PDF model to (x, y) data.

    Parameters
    ----------
    x, y : array-like
    p0 : (mu0, sigma0), optional
        Initial parameter guess.  Defaults to ``(log(median(x)), 1.0)``.
    **kwargs
        Forwarded to :func:`scipy.optimize.curve_fit`.

    Returns
    -------
    popt : ndarray, shape (2,)
        Optimal (μ, σ).
    pcov : ndarray, shape (2, 2)
        Covariance of popt.
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    finite = np.isfinite(x) & (x > 0) & np.isfinite(y)
    if p0 is None:
        p0 = [np.log(np.nanmedian(x[finite])), 1.0]
    return curve_fit(lognormal_forward, x[finite], y[finite], p0=p0, **kwargs)


# ---------------------------------------------------------------------------
# Gaussian (normal) probability density model
# ---------------------------------------------------------------------------

def gaussian_forward(x, mu, sigma):
    """Gaussian PDF forward model: y = norm.pdf(x, mu, sigma)."""
    return norm.pdf(x, loc=mu, scale=sigma)


def gaussian_jac(x, mu, sigma):
    """Jacobian of gaussian_forward w.r.t. (mu, sigma).

    Returns an (N, 2) array: [dy/dmu, dy/dsigma].
    """
    x = np.asarray(x, dtype=float)
    var = sigma**2
    u = x - mu
    arg = np.exp(-0.5 * u**2 / var) / np.sqrt(2 * np.pi)
    dydmu = u * arg / (sigma * var)
    dydsigma = (u**2 - var) * arg / var**2
    return np.column_stack([dydmu, dydsigma])


def fit_gaussian(x, y, p0=None, **kwargs):
    """Fit a Gaussian PDF model to (x, y) data.

    Parameters
    ----------
    x, y : array-like
    p0 : (mu0, sigma0), optional
        Initial parameter guess. Defaults to the mean/std of `x` weighted
        by `y` (i.e. estimated directly from where the density is
        concentrated, not from `x` alone -- using `mean(x)` as the mu0
        guess fails whenever `x` happens to be symmetric about a value
        other than the true peak).
    **kwargs
        Forwarded to :func:`scipy.optimize.curve_fit`.

    Returns
    -------
    popt : ndarray, shape (2,)
        Optimal (mu, sigma).
    pcov : ndarray, shape (2, 2)
        Covariance of popt.
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    if p0 is None:
        w = np.clip(y, 0, None)
        if w.sum() > 0:
            mu0 = np.average(x, weights=w)
            sigma0 = np.sqrt(np.average((x - mu0)**2, weights=w)) or 1.0
        else:
            mu0, sigma0 = np.nanmean(x), np.nanstd(x) or 1.0
        p0 = [mu0, sigma0]
    return curve_fit(gaussian_forward, x, y, p0=p0, **kwargs)


# ---------------------------------------------------------------------------
# Maxwell-Boltzmann probability density model
# ---------------------------------------------------------------------------

def maxwell_forward(x, a):
    """Maxwell-Boltzmann PDF forward model: y = maxwell.pdf(x, scale=a)."""
    return maxwell.pdf(x, scale=a)


def maxwell_jac(x, a):
    """Jacobian of maxwell_forward w.r.t. the scale parameter a.

    Returns an (N, 1) array: [dy/da].
    """
    x = np.asarray(x, dtype=float)
    dyda = np.sqrt(2 / np.pi) * a**-6 * x**2 * (x**2 - 3 * a**2) * np.exp(-0.5 * x**2 / a**2)
    return dyda.reshape(-1, 1)


def fit_maxwell(x, y, p0=None, **kwargs):
    """Fit a Maxwell-Boltzmann PDF model to (x, y) data.

    Parameters
    ----------
    x, y : array-like
    p0 : (a0,), optional
        Initial scale-parameter guess. Defaults to an estimate from the
        mean of `x` weighted by `y` (the density's concentration), using
        the Maxwell-Boltzmann mean-to-scale relation ``mean = 2*a*sqrt(2/pi)``.
    **kwargs
        Forwarded to :func:`scipy.optimize.curve_fit`.

    Returns
    -------
    popt : ndarray, shape (1,)
        Optimal scale parameter a.
    pcov : ndarray, shape (1, 1)
        Covariance of popt.
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    if p0 is None:
        w = np.clip(y, 0, None)
        if w.sum() > 0:
            mean0 = np.average(x, weights=w)
            a0 = mean0 / (2 * np.sqrt(2 / np.pi))
        else:
            a0 = np.nanstd(x) or 1.0
        p0 = [a0 or 1.0]
    return curve_fit(maxwell_forward, x, y, p0=p0, **kwargs)


# ---------------------------------------------------------------------------
# Thermal conductivity vs. anorthite content and temperature
# ---------------------------------------------------------------------------

def tc_forward(X, m1, m2, m3, m4):
    """Thermal conductivity forward model.

    k = (m1 + m2*An + m3*An^2) * (298 / T)^m4

    Parameters
    ----------
    X : array-like, shape (2, N)
        Stacked independent variables ``[An, T]``: molar fraction
        anorthite content and temperature in Kelvin.
    m1, m2, m3 : float
        Quadratic-in-An coefficients for thermal conductivity at 298 K.
    m4 : float
        Temperature power-law exponent.
    """
    an, T = X
    return (m1 + m2 * an + m3 * an**2) * (298.0 / T)**m4


def tc_jac(X, m1, m2, m3, m4):
    """Jacobian of tc_forward w.r.t. (m1, m2, m3, m4).

    Returns an (N, 4) array.
    """
    an, T = X
    base = (298.0 / T)**m4
    dm1 = base
    dm2 = an * base
    dm3 = an**2 * base
    dm4 = (m1 + m2 * an + m3 * an**2) * base * np.log(298.0 / T)
    return np.column_stack([dm1, dm2, dm3, dm4])


def fit_tc(an, T, k, p0=None, **kwargs):
    """Fit the thermal conductivity model to (An, T, k) data.

    Parameters
    ----------
    an : array-like
        Molar fraction anorthite content.
    T : array-like
        Temperature, K.
    k : array-like
        Thermal conductivity, W m^-1 K^-1.
    p0 : (m1, m2, m3, m4), optional
        Initial parameter guess. Defaults to ``(mean(k), 0, 0, 1)``.
    **kwargs
        Forwarded to :func:`scipy.optimize.curve_fit`.

    Returns
    -------
    popt : ndarray, shape (4,)
        Optimal (m1, m2, m3, m4).
    pcov : ndarray, shape (4, 4)
        Covariance of popt.
    """
    an, T, k = np.asarray(an, float), np.asarray(T, float), np.asarray(k, float)
    X = np.vstack([an, T])
    if p0 is None:
        p0 = [np.nanmean(k), 0.0, 0.0, 1.0]
    return curve_fit(tc_forward, X, k, p0=p0, **kwargs)


def tc_plot(an, T, k, popt=None, axes=None):
    """Plots thermal conductivity vs. anorthite content and temperature.

    Left panel shows k vs. An at 298 K; right panel shows k vs. T colored
    by An. If `popt` (from :func:`fit_tc`) is given, the fitted model curve
    is overlaid on the left panel at 298 K.

    Parameters
    ----------
    an : array-like
        Molar fraction anorthite content.
    T : array-like
        Temperature, K.
    k : array-like
        Thermal conductivity, W m^-1 K^-1.
    popt : array-like, shape (4,), optional
        Fitted (m1, m2, m3, m4) from :func:`fit_tc`, overlaid at T=298 K.
    axes : (matplotlib.axes.Axes, matplotlib.axes.Axes), optional
        Existing (left, right) axes to draw into; a new figure is created
        if not given.

    Returns
    -------
    (matplotlib.axes.Axes, matplotlib.axes.Axes)
    """
    import matplotlib.pyplot as plt

    an, T, k = np.asarray(an, float), np.asarray(T, float), np.asarray(k, float)

    if axes is None:
        fig, axes = plt.subplots(1, 2, figsize=(9, 4.5))
    ax1, ax2 = axes

    ind = T == 298
    sc1 = ax1.scatter(an[ind], k[ind], c=an[ind], s=24, cmap='winter')
    if popt is not None:
        an_line = np.linspace(0, 1, 100)
        k_line = tc_forward(np.vstack([an_line, np.full_like(an_line, 298.0)]), *popt)
        ax1.plot(an_line, k_line, 'k-')
    plt.colorbar(sc1, ax=ax1, label='An')
    ax1.set_xlabel('% Anorthite')
    ax1.set_ylabel('Thermal Conductivity [W/m/K]')
    ax1.set_xlim(0, 1)
    ax1.set_ylim(1, 3)
    ax1.set_box_aspect(1)

    sc2 = ax2.scatter(T, k, c=an, s=24, cmap='winter')
    plt.colorbar(sc2, ax=ax2, label='An')
    ax2.set_xlabel('Temperature [K]')
    ax2.set_ylabel('Thermal conductivity [W/m/K]')
    ax2.set_xlim(250, 700)
    ax2.set_ylim(1, 3)
    ax2.set_box_aspect(1)

    return ax1, ax2
