import numpy as np
from scipy.optimize import brentq

"""Useful tools for kernel density estimation (KDE) and its derivatives, including
peak location and uncertainty estimation. The main function is peak_uncertainty(), which returns a dict with the peak location, its standard error, a 95% CI, and the fraction of bootstrap samples that produced a valid peak. The KDE is Gaussian, with bandwidth h.

Example:
    m   = find_mode(data, h, lo, hi)
    f0  = kde_deriv(m, data, h,          0)[0]
    f2  = kde_deriv(m, data, 1.5*h,      2)[0]     # pilot bandwidth for curvature
    se  = np.sqrt(f0 / (4*np.sqrt(np.pi) * data.size * h**3)) / abs(f2)
"""

def kde_deriv(x0, data, h, order=0):
    """Gaussian KDE and its derivatives, evaluated at scalar/array x0."""
    z   = (np.atleast_1d(x0)[:, None] - data[None, :]) / h
    phi = np.exp(-0.5*z**2) / np.sqrt(2*np.pi)
    n   = data.size
    if order == 0: return  phi.sum(1) / (n*h)
    if order == 1: return -(z*phi).sum(1) / (n*h**2)
    if order == 2: return ((z**2 - 1)*phi).sum(1) / (n*h**3)

def find_mode(data, h, lo, hi, ngrid=512):
    """Locate the peak in [lo, hi]: coarse grid, then root-find on f'."""
    g  = np.linspace(lo, hi, ngrid)
    f  = kde_deriv(g, data, h, 0)
    i  = np.argmax(f)
    if i in (0, ngrid-1):
        return np.nan                      # peak ran off the bracket
    return brentq(lambda t: kde_deriv(t, data, h, 1)[0], g[i-1], g[i+1])

def peak_uncertainty(data, h, lo, hi, B=2000, sigma=None, smoothed=False, seed=0):
    """Bootstrap CI for a KDE peak position. sigma: per-point 1s analytical errors."""
    rng   = np.random.default_rng(seed)
    n     = data.size
    peaks = np.empty(B)
    for b in range(B):
        idx = rng.integers(0, n, n)        # resample -> sampling uncertainty
        xb  = data[idx]
        if sigma is not None:              # perturb -> measurement uncertainty
            xb = rng.normal(xb, sigma[idx])
        if smoothed:                       # smoothed bootstrap, variance-corrected
            s2 = xb.var(ddof=1)
            xb = xb.mean() + (xb - xb.mean() + h*rng.standard_normal(n)) \
                 / np.sqrt(1 + h**2/s2)
        peaks[b] = find_mode(xb, h, lo, hi)

    ok = np.isfinite(peaks)
    return dict(peak       = find_mode(data, h, lo, hi),
                se         = np.nanstd(peaks, ddof=1),
                ci95       = np.nanpercentile(peaks[ok], [2.5, 97.5]),
                detect_frac= ok.mean())    # how often the peak exists at all