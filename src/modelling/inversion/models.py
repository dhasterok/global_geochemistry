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
