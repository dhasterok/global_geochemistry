"""
Radiogenic heat production estimation from whole-rock geochemistry.

Two computation methods are available:

``'mcdonough'`` (default)
    McDonough et al. (2020), Geochem. Geophys. Geosyst. 21, e2019GC008865.
    Equation given in that paper (all concentrations in ppm)::

        h [μW kg⁻¹] = 10⁻⁶ × (98.29 C_U + 26.18 C_Th
                                + 3.387×10⁻³ C_K
                                + 0.01139 C_Rb + 0.04595 C_Sm)

    Rb and Sm are optional (combined ~0.2% of UCC heat production).

``'rybach'``
    Rybach (1988) in Haenel et al. (eds), Handbook of Terrestrial Heat-Flow.
    ::

        H [μW kg⁻¹] = 10⁻⁵ × (9.52 C_U + 2.56 C_Th + 3.48 C_K)

    where C_U, C_Th in ppm and C_K in wt% elemental K (not K₂O).

BDL convention in input DataFrame columns
------------------------------------------
    ``x > 0``   measured above detection limit — used directly.
    ``x < 0``   below detection limit (BDL); ``|x|`` is the detection limit.
    ``x = 0``   BDL with unknown detection limit.
    ``x = NaN`` not measured; optionally predicted from correlated elements.

Missing-value prediction (``estimate_missing=True``)
----------------------------------------------------
    U ↔ Th    via rock-type median log₁₀(Th/U).
    Rb        from K (K/Rb ≈ 300 for continental crust).
    Sm        from Nd (Sm/Nd ≈ 0.22 for upper continental crust).

Gaussian censoring (``censor_bdl=True``)
----------------------------------------
    Uses the Michael & Schucany (1986) method implemented in
    ``src.processing.censored.gausscensor`` to fit a log-normal distribution
    to each element, then substitutes the conditional expectation
    E[log₁₀ X | X ≤ d] for each BDL observation with detection limit *d*.

Interval censoring
------------------
    Trace-element concentrations are often reported to 2–3 significant figures,
    creating apparent interval censoring.  With more than ~10 values per decade
    the effect on heat-production statistics is negligible compared to typical
    5–10% analytical uncertainty and is not corrected here.

References
----------
McDonough, W.F., Šrámek, O., & Mound, J.E. (2020) Radiogenic heating and
    its influence on Rocky Planet dynamics. Geochem. Geophys. Geosyst. 21,
    e2019GC008865. https://doi.org/10.1029/2019GC008865

Rybach, L. (1988) Determination of heat production rate. In: Haenel R.,
    Rybach L., Stegena L. (eds) Handbook of Terrestrial Heat-Flow Density
    Determination. Kluwer, Dordrecht. pp. 125–142.
"""

import warnings
import numpy as np
import pandas as pd
from scipy.stats import norm as _norm
from global_geochemistry.utils.molecular import MolecularWeightCalculator

_mwc = MolecularWeightCalculator()

# ---------------------------------------------------------------------------
# Molecular-weight constants
# ---------------------------------------------------------------------------

_K2O_TO_K = 2 * _mwc.molecular_weight('K') / _mwc.molecular_weight('K2O')

# ---------------------------------------------------------------------------
# Heat-production coefficients
# ---------------------------------------------------------------------------
# Units stored for each element:
#   u, th, rb, sm : μW kg⁻¹ per ppm of the natural element
#   k              : μW kg⁻¹ per wt% of elemental K  (NOT K₂O)
#
# McDonough (2020) formula in ppm for all elements:
#   h = 1e-6 × (98.29 U + 26.18 Th + 3.387e-3 K_ppm + 0.01139 Rb + 0.04595 Sm)
# Converted to per wt% K: 3.387e-3 × 1e-6 × 1e4 wt%⁻¹ = 3.387e-5

_HP_CONSTANTS = {
    'rybach': {
        'u' : 9.52e-5,
        'th': 2.56e-5,
        'k' : 3.48e-5,
        'rb': 0.0,
        'sm': 0.0,
        '_ref': 'Rybach (1988)',
    },
    'mcdonough': {
        'u' : 9.829e-5,
        'th': 2.618e-5,
        'k' : 3.387e-5,
        'rb': 1.139e-8,
        'sm': 4.595e-8,
        '_ref': 'McDonough et al. (2020) GC008865',
    },
}

# Legacy Rybach constants kept for compute_hp() backward compatibility
_CK    = 3.48
_CU    = 9.52
_CTH   = 2.56
_SCALE = 1e-5


# ---------------------------------------------------------------------------
# Time-dependent heat production — nuclear constants
# ---------------------------------------------------------------------------

# Age of Solar System formation (Bouvier & Wadhwa 2010)
_T_SS = 4568.0  # Ma

# Per-isotope data for time-correction of long-lived nuclides (McDonough 2020).
# Tuple: (half_life_Ma, H_per_kg_isotope [W/kg], present-day mass-fraction in element)
# The mass fractions are the PRESENT-DAY (t=0) isotopic abundances.
# K, Rb, Sm: only one isotope contributes meaningfully to heat.
# U: two isotopes (²³⁸U and ²³⁵U) with very different half-lives — must be kept separate.
_ISO_LONG = {
    '238U':  (4468.0,    9.474e-5,  0.992742),  # fraction of natural U
    '235U':  (703.8,     5.688e-4,  0.007204),  # fraction of natural U
    '232Th': (14050.0,   2.638e-5,  1.000000),  # essentially monoisotopic
    '40K':   (1248.0,    2.902e-5,  1.167e-4),  # fraction of natural K
    '87Rb':  (49700.0,   4.09e-8,   0.27835),   # fraction of natural Rb
    '147Sm': (106000.0,  3.063e-7,  0.15000),   # fraction of natural Sm
}

# Short-lived nuclides for Hadean / early Solar System calculations.
# These are fully extinct in the modern Earth.  Their initial concentrations
# are inferred from canonical abundance ratios at Solar System formation (t_SS).
#
# Tuple: (half_life_Ma, H_per_kg_isotope [W/kg], default_initial_ratio,
#         mass_ratio_to_stable, f_stable_in_element, source_element_column)
#
# Initial ratio definitions:
#   26Al : (²⁶Al/²⁷Al)_0 × (26/27)  → ppm ²⁶Al from al_ppm
#   60Fe : (⁶⁰Fe/⁵⁶Fe)_0 × (60/56) × f_56Fe → ppm ⁶⁰Fe from fe_ppm
#   10Be : (¹⁰Be/⁹Be)_0  × (10/9)   → ppm ¹⁰Be from be_ppm
#   244Pu: (²⁴⁴Pu/²³⁸U)_0 × (244/238) × f_238U → ppm ²⁴⁴Pu from u_ppm
#   146Sm: (¹⁴⁶Sm/¹⁴⁴Sm)_0 × (146/144) × f_144Sm → ppm ¹⁴⁶Sm from sm_ppm
#
# H values computed from λ × N_A/M × Q_deposited.  Q values from literature:
#   26Al: ~3.12 MeV deposited (Castillo-Rogez et al. 2009)
#   60Fe: ~3.56 MeV (two-step β chain through ⁶⁰Co; Moskovitz & Gaidos 2011)
#   10Be: ~0.556 MeV β (most deposited)
#   244Pu: ~5.80 MeV α decay (dominant channel)
#   146Sm: ~2.53 MeV α decay
#
# Initial ratios (default): canonical Solar System values.
# ⁶⁰Fe/⁵⁶Fe_0 = 1e-6 is the upper estimate; recent meteoritic data range from
# 1e-8 to 1e-6 — this is the most uncertain parameter here.
_ISO_SHORT = {
    '26Al':  (0.717,   0.355,   5.23e-5, 26/27,     1.0000, 'al_ppm'),
    '60Fe':  (2.62,    0.048,   1.0e-6,  60/56,     0.9172, 'fe_ppm'),
    '10Be':  (1.387,   0.028,   7.0e-4,  10/9,      1.0000, 'be_ppm'),
    '244Pu': (81.3,    6.20e-4, 6.8e-3,  244/238,   0.9927, 'u_ppm'),
    '146Sm': (103.0,   3.56e-4, 8.2e-3,  146/144,   0.0308, 'sm_ppm'),
}


def _hp_time_factors(t):
    """Element-level time-correction factors at time *t* (Ma, positive = past).

    Returns a dict with keys ``'u'``, ``'th'``, ``'k'``, ``'rb'``, ``'sm'``.
    Multiplying the present-day coefficient for each element by its factor gives
    the coefficient that should be used at time *t*.

    The U factor accounts for the different decay rates of ²³⁸U and ²³⁵U using
    their present-day isotopic abundances and per-isotope heat production rates.
    At t = 0 every factor equals 1.

    Parameters
    ----------
    t : float or array-like
        Time in Ma.  Positive = past; negative = future.

    Returns
    -------
    dict of numpy.ndarray
    """
    t = np.asarray(t, dtype=float)

    # U: two components with very different half-lives
    t12_238, H_238, f_238 = _ISO_LONG['238U']
    t12_235, H_235, f_235 = _ISO_LONG['235U']
    l238 = np.log(2) / t12_238   # 1/Ma
    l235 = np.log(2) / t12_235

    c238_0 = f_238 * H_238
    c235_0 = f_235 * H_235
    c_U0   = c238_0 + c235_0     # denominator (= present-day combined)

    gU = (c238_0 * np.exp(l238 * t) + c235_0 * np.exp(l235 * t)) / c_U0

    # Single-isotope elements
    gTh = np.exp(np.log(2) / _ISO_LONG['232Th'][0] * t)
    gK  = np.exp(np.log(2) / _ISO_LONG['40K'][0]   * t)
    gRb = np.exp(np.log(2) / _ISO_LONG['87Rb'][0]  * t)
    gSm = np.exp(np.log(2) / _ISO_LONG['147Sm'][0] * t)

    return {'u': gU, 'th': gTh, 'k': gK, 'rb': gRb, 'sm': gSm}


def _short_lived_hp(t, t_ss=_T_SS, *,
                    Al=None, Fe=None, Be=None,
                    Sm=None, U=None,
                    al26_al27_0=None, fe60_fe56_0=None,
                    be10_be9_0=None, sm146_sm144_0=None,
                    pu244_u238_0=None):
    """Heat production from short-lived nuclides at time *t* (μW kg⁻¹).

    Concentrations at time *t* are computed as::

        C_i(t) = C_i(t_ss) × exp(−λ_i × (t_ss − t))

    where ``C_i(t_ss)`` is inferred from the canonical initial abundance ratio
    (see ``_ISO_SHORT``).  Returns zero for t ≤ 0 (present/future, where
    all short-lived nuclides are effectively extinct).

    Parameters
    ----------
    t : float or array-like
        Age in Ma (positive = past).  Scalar or per-sample array.
    t_ss : float
        Solar System formation age in Ma.  Default 4568.0.
    Al, Fe, Be : array-like or None
        Present-day Al, Fe, Be concentrations in ppm.  Required for ²⁶Al,
        ⁶⁰Fe, ¹⁰Be contributions respectively.
    Sm, U : array-like or None
        Present-day Sm and U in ppm, required for ¹⁴⁶Sm and ²⁴⁴Pu.
    al26_al27_0, fe60_fe56_0, be10_be9_0, sm146_sm144_0, pu244_u238_0 : float or None
        Override canonical initial abundance ratios.

    Returns
    -------
    numpy.ndarray
        Short-lived heat production in μW kg⁻¹.
    """
    t = np.asarray(t, dtype=float)

    # Determine output shape from the first non-None element array
    ref = next((a for a in (Al, Fe, Be, Sm, U) if a is not None), None)
    if ref is None:
        return np.zeros_like(t)
    h = np.zeros_like(np.asarray(ref, dtype=float))

    def _contrib(iso_key, C, ratio_override=None):
        """Add contribution from one short-lived nuclide.

        C : array-like, present-day element concentration in ppm.
        H [W/kg_isotope] × C [ppm] = μW/kg_rock  (1e-6 × 1e6 cancel).
        """
        t12, H, ratio0, mass_ratio, f_stable, _ = _ISO_SHORT[iso_key]
        lam = np.log(2) / t12   # 1/Ma
        ratio = ratio_override if ratio_override is not None else ratio0

        # Initial concentration of short-lived isotope at t_ss [ppm]
        C0 = ratio * mass_ratio * f_stable * np.asarray(C, dtype=float)

        # Decay from t_ss to t:  exp(−λ × (t_ss − t))
        # t ≤ 0  → contribution = 0 (all extinct by present day)
        # t > t_ss → clamp to t_ss (concentration at formation)
        dt = np.clip(t_ss - np.asarray(t, dtype=float), 0.0, None)
        decay = np.where(t > 0, np.exp(-lam * dt), 0.0)

        # Broadcast C0 × decay (C0 may be (N,), decay may be scalar or (N,))
        return H * C0 * decay

    if Al is not None:
        h = h + _contrib('26Al',  Al, al26_al27_0)
    if Fe is not None:
        h = h + _contrib('60Fe',  Fe, fe60_fe56_0)
    if Be is not None:
        h = h + _contrib('10Be',  Be, be10_be9_0)
    if U is not None:
        h = h + _contrib('244Pu', U,  pu244_u238_0)
    if Sm is not None:
        h = h + _contrib('146Sm', Sm, sm146_sm144_0)

    return h


# ---------------------------------------------------------------------------
# Low-level heat-production computation
# ---------------------------------------------------------------------------

def radiogenic_hp(U, Th, K, Rb=None, Sm=None, method='mcdonough'):
    """Compute specific heat production (μW kg⁻¹) from element concentrations.

    All inputs must be *positive* (detected) values; apply BDL handling and
    missing-value prediction before calling this function.  NaN propagates
    through the calculation.

    Parameters
    ----------
    U, Th : array-like
        Uranium and thorium concentrations in ppm.
    K : array-like
        Elemental potassium in wt% (**not** K₂O).
    Rb, Sm : array-like or None
        Rubidium and samarium in ppm.  Ignored unless *method* is ``'mcdonough'``.
    method : {'mcdonough', 'rybach'}

    Returns
    -------
    numpy.ndarray
        Heat production in μW kg⁻¹.
    """
    if method not in _HP_CONSTANTS:
        raise ValueError(f"method must be one of {list(_HP_CONSTANTS)}; got {method!r}")
    c = _HP_CONSTANTS[method]

    U  = np.asarray(U,  float)
    Th = np.asarray(Th, float)
    K  = np.asarray(K,  float)

    h = c['u'] * U + c['th'] * Th + c['k'] * K

    if Rb is not None and c['rb'] > 0:
        h = h + c['rb'] * np.asarray(Rb, float)
    if Sm is not None and c['sm'] > 0:
        h = h + c['sm'] * np.asarray(Sm, float)

    return h


def radiogenic_hp_at_time(t, U, Th, K,
                           Rb=None, Sm=None,
                           method='mcdonough',
                           early_earth=False,
                           Al=None, Fe=None, Be=None,
                           al26_al27_0=None, fe60_fe56_0=None,
                           be10_be9_0=None, sm146_sm144_0=None,
                           pu244_u238_0=None,
                           t_ss=_T_SS):
    """Compute specific heat production (μW kg⁻¹) at time *t*.

    Each long-lived isotope's contribution is scaled by ``exp(λᵢ × t)``
    relative to its present-day value, using present-day isotopic abundances
    and per-isotope heat-production rates.  The U correction accounts
    separately for ²³⁸U and ²³⁵U.

    Optionally includes short-lived extinct nuclides that matter in the
    Hadean and early Solar System (see ``early_earth``).

    Parameters
    ----------
    t : float or array-like
        Time in Ma.  **Positive = past, negative = future.**
        May be a scalar (same time for all samples) or a per-sample array
        (e.g. ``data['avg_age'].to_numpy()``).
    U, Th : array-like
        Present-day uranium and thorium concentrations in ppm.
    K : array-like
        Present-day elemental potassium in wt% (**not** K₂O).
    Rb, Sm : array-like or None
        Present-day rubidium and samarium in ppm.  Only used for
        ``method='mcdonough'``.
    method : {'mcdonough', 'rybach'}
    early_earth : bool
        Include short-lived nuclides: ²⁶Al, ⁶⁰Fe, ¹⁰Be, ²⁴⁴Pu, ¹⁴⁶Sm.
        Only significant for t ≳ 3900 Ma (¹⁴⁶Sm, ²⁴⁴Pu) and
        t ≳ 4500 Ma (²⁶Al, ⁶⁰Fe, ¹⁰Be).  Al, Fe, Be must be supplied.
    Al : array-like or None
        Aluminium in ppm — required for ²⁶Al.
    Fe : array-like or None
        Iron in ppm — required for ⁶⁰Fe.
    Be : array-like or None
        Beryllium in ppm — required for ¹⁰Be.
    al26_al27_0 : float or None
        Override canonical ²⁶Al/²⁷Al initial ratio (default 5.23 × 10⁻⁵).
    fe60_fe56_0 : float or None
        Override canonical ⁶⁰Fe/⁵⁶Fe initial ratio (default 1.0 × 10⁻⁶;
        highly uncertain — published estimates span 10⁻⁸ to 10⁻⁶).
    be10_be9_0 : float or None
        Override canonical ¹⁰Be/⁹Be initial ratio (default 7.0 × 10⁻⁴).
    sm146_sm144_0 : float or None
        Override canonical ¹⁴⁶Sm/¹⁴⁴Sm initial ratio (default 8.2 × 10⁻³).
    pu244_u238_0 : float or None
        Override canonical ²⁴⁴Pu/²³⁸U initial ratio (default 6.8 × 10⁻³).
    t_ss : float
        Solar System formation age in Ma.  Default 4568.0.

    Returns
    -------
    numpy.ndarray
        Heat production in μW kg⁻¹.

    Notes
    -----
    Time-correction factors are defined so that the present-day (t = 0)
    result is identical to ``radiogenic_hp``.  For future times (t < 0)
    all factors are < 1.

    Short-lived nuclide H values (W kg⁻¹ of isotope) are calculated from
    λ N_A/M × Q_deposited using the energies listed in ``_ISO_SHORT``; there
    is ~15% uncertainty in Q for ²⁶Al depending on what fraction of the γ
    ray energy is retained.  For ⁶⁰Fe the initial ⁶⁰Fe/⁵⁶Fe ratio is the
    dominant uncertainty (two orders of magnitude spread in the literature).
    """
    if method not in _HP_CONSTANTS:
        raise ValueError(f"method must be one of {list(_HP_CONSTANTS)}; got {method!r}")

    g = _hp_time_factors(t)
    c = _HP_CONSTANTS[method]

    U  = np.asarray(U,  dtype=float)
    Th = np.asarray(Th, dtype=float)
    K  = np.asarray(K,  dtype=float)

    h = c['u'] * g['u'] * U + c['th'] * g['th'] * Th + c['k'] * g['k'] * K

    if Rb is not None and c['rb'] > 0:
        h = h + c['rb'] * g['rb'] * np.asarray(Rb, dtype=float)
    if Sm is not None and c['sm'] > 0:
        h = h + c['sm'] * g['sm'] * np.asarray(Sm, dtype=float)

    if early_earth:
        h = h + _short_lived_hp(
            t, t_ss=t_ss,
            Al=Al, Fe=Fe, Be=Be,
            Sm=Sm, U=U,
            al26_al27_0=al26_al27_0, fe60_fe56_0=fe60_fe56_0,
            be10_be9_0=be10_be9_0, sm146_sm144_0=sm146_sm144_0,
            pu244_u238_0=pu244_u238_0,
        )

    return h


def compute_hp(k2o, u, th):
    """Compute specific heat production (μW kg⁻¹) — Rybach (1988).

    Parameters
    ----------
    k2o : array-like
        K₂O concentration (wt%).
    u, th : array-like
        U and Th concentrations (ppm).

    Returns
    -------
    numpy.ndarray
        Heat production in μW kg⁻¹.
    """
    k2o = np.asarray(k2o, float)
    u   = np.asarray(u,   float)
    th  = np.asarray(th,  float)
    k   = _K2O_TO_K * np.abs(k2o)
    return _SCALE * (_CK * k + _CU * np.abs(u) + _CTH * np.abs(th))


# ---------------------------------------------------------------------------
# BDL censoring helper
# ---------------------------------------------------------------------------

def _bdl_censor(x, verbose=False, label='element'):
    """Replace BDL values with conditional log-normal expectations.

    Uses ``gausscensor`` to fit a log-normal distribution to the full
    array (detected + censored), then substitutes E[log₁₀ X | X ≤ d]
    (truncated-normal conditional mean) for each BDL sample.

    Parameters
    ----------
    x : numpy.ndarray
        Measured values following the BDL flag convention:
        positive = detected, negative = BDL with |x| as detection limit,
        zero = BDL with unknown detection limit.

    Returns
    -------
    numpy.ndarray
        Copy of *x* with BDL values replaced by positive estimates.
    """
    x = np.asarray(x, dtype=float).copy()
    is_bdl = np.isfinite(x) & (x <= 0)
    if not is_bdl.any():
        return x

    # Unknown DL (zero): substitute median of known DLs
    known_dls = np.abs(x[x < 0])
    default_dl = (np.median(known_dls) if len(known_dls) > 0
                  else np.nanmin(x[x > 0]) / 10 if (x > 0).any()
                  else 0.01)
    zero_mask = np.isfinite(x) & (x == 0)
    if zero_mask.any():
        x[zero_mask] = -default_dl

    # Fit log-normal using gausscensor
    mu = sigma = np.nan
    try:
        from global_geochemistry.geochem.processing.censored import gausscensor
        result = gausscensor(x, scale='log10')
        if isinstance(result, tuple):
            stats = result[0]
        else:
            stats = result
        mu    = float(stats['mu'])
        sigma = float(stats['sigma'])
    except Exception:
        pass

    x_out = x.copy()
    for i in np.where(is_bdl)[0]:
        d = abs(x[i])
        if np.isnan(mu) or np.isnan(sigma) or sigma <= 0 or d <= 0:
            x_out[i] = d / np.sqrt(2)
            continue

        log_d = np.log10(d)
        z_d   = (log_d - mu) / sigma
        phi_z = _norm.pdf(z_d)
        Phi_z = _norm.cdf(z_d)

        if Phi_z > 1e-10:
            # Conditional expectation E[log10 X | log10 X ≤ log10 d]
            E_log_x = mu - sigma * phi_z / Phi_z
            x_out[i] = 10.0 ** E_log_x
        else:
            x_out[i] = d / np.sqrt(2)

    if verbose:
        n = is_bdl.sum()
        print(f'  hpest: {n} BDL {label} values replaced by censored estimates')

    return x_out


# ---------------------------------------------------------------------------
# Element extraction helpers
# ---------------------------------------------------------------------------

def _col_or_nan(data, name):
    """Return column as float array, or NaN array if absent."""
    if name in data.columns:
        return data[name].to_numpy(float)
    return np.full(len(data), np.nan)


def get_k2o(data):
    """Return K₂O (wt%) as a numpy array, converting from K_ppm if needed."""
    if 'k2o' in data.columns:
        return data['k2o'].to_numpy(float)
    if 'k_ppm' in data.columns:
        factor = _mwc.molecular_weight('K2O') / (2 * _mwc.molecular_weight('K') * 1e4)
        return data['k_ppm'].to_numpy(float) * factor
    return np.full(len(data), np.nan)


# ---------------------------------------------------------------------------
# Main HP estimator
# ---------------------------------------------------------------------------

def hpest(data, method='mcdonough', include_rb=False, include_sm=False,
          estimate_missing=True, censor_bdl=True,
          bdl_k2o=-0.01, bdl_u=-1.0, bdl_th=-1.0,
          t=None, early_earth=False,
          verbose=True):
    """Estimate radiogenic heat production from whole-rock geochemistry.

    Adds ``heat_production_mass`` (μW kg⁻¹) and ``heat_production``
    (μW m⁻³) columns.  Density is taken from ``density_model`` if present.

    Parameters
    ----------
    data : pandas.DataFrame
        Must contain ``u_ppm``, ``th_ppm``, and either ``k2o`` or ``k_ppm``.
        Optionally ``rb_ppm``, ``sm_ppm``, ``nd_ppm``.
    method : {'mcdonough', 'rybach'}
        Heat-production formula.  ``'mcdonough'`` (default) uses updated
        decay constants from McDonough et al. (2020).
    include_rb : bool
        Include Rb contribution (McDonough only, default False).
    include_sm : bool
        Include Sm contribution (McDonough only, default False).
    estimate_missing : bool
        If True (default), predict missing or BDL concentrations from
        correlated elements before computing HP:

        - U ↔ Th  via rock-type median log₁₀(Th/U).
        - Rb from K (K/Rb ≈ 300 for continental crust).
        - Sm from Nd (Sm/Nd ≈ 0.22 for upper continental crust).

        Predictions are only applied where the target element is missing/BDL
        *and* the source element is positively measured.
    censor_bdl : bool
        If True (default), remaining BDL values after prediction are
        substituted with a Gaussian-censoring estimate (conditional
        log-normal expectation via ``gausscensor``).  If False, the
        absolute value of the detection limit is used directly.
    bdl_k2o, bdl_u, bdl_th : float
        Substitute detection limit for elements whose reported value is
        exactly zero (unknown DL).  Negative values flag BDL.
    t : float, array-like, or None
        Time in Ma at which to compute HP (positive = past, negative = future).
        If None (default), computes present-day HP using ``radiogenic_hp``.
        May be a scalar (single time for all samples) or a per-sample array
        the same length as *data* (e.g. ``data['avg_age'].to_numpy()``).
        When set, uses isotope-decay-corrected coefficients via
        ``radiogenic_hp_at_time``.
    early_earth : bool
        Include short-lived extinct nuclides (²⁶Al, ⁶⁰Fe, ¹⁰Be, ²⁴⁴Pu,
        ¹⁴⁶Sm).  Only meaningful for t ≳ 3900 Ma (¹⁴⁶Sm, ²⁴⁴Pu) and
        t ≳ 4400 Ma (²⁶Al, ⁶⁰Fe, ¹⁰Be).  Requires ``t`` to be set.
        Al and Fe are extracted from ``al2o3`` and ``feo_tot`` columns and
        converted to elemental ppm; Be is taken from ``be_ppm`` if present.
    verbose : bool

    Returns
    -------
    pandas.DataFrame
        Modified copy of *data* with columns added:

        heat_production_mass   μW kg⁻¹ (mass-specific)
        heat_production        μW m⁻³  (volumetric; requires density_model)
        heat_production_bdl    bool — True where any BDL substitution was made
    """
    if method not in _HP_CONSTANTS:
        raise ValueError(f"method must be one of {list(_HP_CONSTANTS)}; got {method!r}")

    data = data.copy()
    n    = len(data)

    if verbose:
        c = _HP_CONSTANTS[method]
        print(f'  hpest: method = {method!r}  ({c["_ref"]})')

    # ------------------------------------------------------------------ #
    # 1. Extract raw element concentrations                               #
    # ------------------------------------------------------------------ #

    k2o_raw  = get_k2o(data)
    u_raw    = _col_or_nan(data, 'u_ppm')
    th_raw   = _col_or_nan(data, 'th_ppm')
    rb_raw   = _col_or_nan(data, 'rb_ppm') if include_rb else np.full(n, np.nan)
    sm_raw   = _col_or_nan(data, 'sm_ppm') if include_sm else np.full(n, np.nan)
    nd_raw   = _col_or_nan(data, 'nd_ppm')

    # Rock-type masks for ratio-based prediction
    grp = (data['rock_group'].str.lower()
           if 'rock_group' in data.columns else pd.Series([''] * n))
    org = (data['rock_origin'].str.lower()
           if 'rock_origin' in data.columns else pd.Series([''] * n))
    ign_mask = ((grp == 'igneous') |
                org.isin(['metaigneous', 'metavolcanic', 'metaplutonic'])).to_numpy()
    sed_mask = ((grp == 'sedimentary') |
                (org == 'metasedimentary')).to_numpy()

    # ------------------------------------------------------------------ #
    # 2. Substitute zero-BDL with default detection limits               #
    # ------------------------------------------------------------------ #

    k2o_raw = np.where(np.isfinite(k2o_raw) & (k2o_raw == 0), bdl_k2o, k2o_raw)
    u_raw   = np.where(np.isfinite(u_raw)   & (u_raw   == 0), bdl_u,   u_raw)
    th_raw  = np.where(np.isfinite(th_raw)  & (th_raw  == 0), bdl_th,  th_raw)

    # Convert K2O → elemental K in wt%, preserving BDL sign
    k_raw = np.where(np.isfinite(k2o_raw),
                     np.sign(k2o_raw) * _K2O_TO_K * np.abs(k2o_raw),
                     np.nan)

    # Working copies (will be modified by prediction / censoring steps)
    u  = u_raw.copy()
    th = th_raw.copy()
    k  = k_raw.copy()
    rb = rb_raw.copy()
    sm = sm_raw.copy()

    bdl_flag = np.zeros(n, dtype=bool)  # track samples that used BDL substitution

    # ------------------------------------------------------------------ #
    # 3. Predict missing / BDL values from correlated elements           #
    # ------------------------------------------------------------------ #

    if estimate_missing:
        _predict_u_th(u, th, u_raw, th_raw, ign_mask, sed_mask, verbose)

        if include_rb:
            _predict_rb_from_k(rb, k, verbose)

        if include_sm:
            _predict_sm_from_nd(sm, nd_raw, verbose)

    # ------------------------------------------------------------------ #
    # 4. Handle remaining BDL values                                      #
    # ------------------------------------------------------------------ #

    for arr, lbl in [(u, 'U'), (th, 'Th'), (k, 'K'), (rb, 'Rb'), (sm, 'Sm')]:
        if arr is None:
            continue
        is_bdl = np.isfinite(arr) & (arr <= 0)
        bdl_flag |= is_bdl
        if not is_bdl.any():
            continue
        if censor_bdl:
            arr[:] = _bdl_censor(arr, verbose=verbose, label=lbl)
        else:
            arr[is_bdl] = np.abs(arr[is_bdl])   # use DL as-is (upper bound)

    # After substitution, any remaining ≤0 become NaN (shouldn't happen)
    for arr in (u, th, k, rb, sm):
        arr[np.isfinite(arr) & (arr <= 0)] = np.nan

    # ------------------------------------------------------------------ #
    # 5. Compute heat production                                          #
    # ------------------------------------------------------------------ #

    rb_in = rb if include_rb else None
    sm_in = sm if include_sm else None

    if t is not None:
        if verbose:
            t_scalar = np.asarray(t)
            if t_scalar.ndim == 0:
                print(f'  hpest: time-corrected HP at {float(t_scalar):.1f} Ma')
            else:
                print(f'  hpest: time-corrected HP (per-sample ages)')

        # Extract major-element concentrations for short-lived nuclides
        Al_in = Fe_in = Be_in = None
        if early_earth:
            if 'al2o3' in data.columns:
                mw_al2o3 = _mwc.molecular_weight('Al2O3')
                mw_al    = _mwc.molecular_weight('Al')
                al2o3    = data['al2o3'].to_numpy(float)
                Al_in    = np.where(np.isfinite(al2o3),
                                    al2o3 * (2 * mw_al / mw_al2o3) * 1e4,
                                    np.nan)
            if 'feo_tot' in data.columns:
                mw_feo = _mwc.molecular_weight('FeO')
                mw_fe  = _mwc.molecular_weight('Fe')
                feo    = data['feo_tot'].to_numpy(float)
                Fe_in  = np.where(np.isfinite(feo),
                                  feo * (mw_fe / mw_feo) * 1e4,
                                  np.nan)
            Be_raw = _col_or_nan(data, 'be_ppm')
            Be_in  = Be_raw if np.isfinite(Be_raw).any() else None

        h_mass = radiogenic_hp_at_time(
            t, u, th, k,
            Rb=rb_in, Sm=sm_in,
            method=method,
            early_earth=early_earth,
            Al=Al_in, Fe=Fe_in, Be=Be_in,
        )
    else:
        h_mass = radiogenic_hp(u, th, k, Rb=rb_in, Sm=sm_in, method=method)

    # NaN where all three core inputs are NaN
    all_nan = np.isnan(u) & np.isnan(th) & np.isnan(k)
    h_mass[all_nan] = np.nan

    data['heat_production_mass'] = h_mass
    data['heat_production_bdl']  = bdl_flag

    # Volumetric HP (μW m⁻³)
    if 'density_model' not in data.columns:
        from global_geochemistry.geochem.physprop.seismic import vpest
        from global_geochemistry.geochem.physprop.density import density_cm
        data = vpest(data)
        data['density_model'] = density_cm(data['p_velocity'])

    data['heat_production'] = h_mass * data['density_model']

    if verbose:
        n_valid = int(np.isfinite(h_mass).sum())
        n_bdl   = int(bdl_flag.sum())
        print(f'  hpest: {n_valid:,} samples with HP; {n_bdl:,} used BDL substitution')

    return data


# ---------------------------------------------------------------------------
# Missing-value prediction helpers (used by hpest)
# ---------------------------------------------------------------------------

_K_RB_RATIO = 300.0   # typical K/Rb for continental crust
_SM_ND_RATIO = 0.22   # Sm/Nd for upper continental crust (Taylor & McLennan 1985)


def _predict_u_th(u, th, u_raw, th_raw, ign_mask, sed_mask, verbose):
    """In-place: predict missing U from Th and vice versa.

    Uses rock-type-specific median log₁₀(Th/U) computed from samples
    that have both values measured (positive).
    """
    sio2 = np.full(len(u), np.nan)  # placeholder; ratio varies by SiO2 for sed

    both = (u_raw > 0) & (th_raw > 0)
    with np.errstate(divide='ignore', invalid='ignore'):
        log_ratio = np.where(both, np.log10(th_raw / u_raw), np.nan)

    def _med(mask):
        r = log_ratio[mask & both]
        r = r[np.isfinite(r)]
        return float(np.median(r)) if len(r) else np.nan

    med_ign = _med(ign_mask)
    med_sed = _med(sed_mask)
    med_all = _med(np.ones(len(u), dtype=bool))  # fallback

    def _ratio(i):
        if ign_mask[i]:
            return med_ign if np.isfinite(med_ign) else med_all
        if sed_mask[i]:
            return med_sed if np.isfinite(med_sed) else med_all
        return med_all

    n_u = n_th = 0
    for i in range(len(u)):
        r = _ratio(i)
        if not np.isfinite(r):
            continue

        # Predict U from Th
        u_missing = not np.isfinite(u[i]) or u[i] <= 0
        th_pos    = np.isfinite(th[i]) and th[i] > 0
        if u_missing and th_pos:
            u[i] = th[i] / 10**r
            n_u += 1

        # Predict Th from U
        th_missing = not np.isfinite(th[i]) or th[i] <= 0
        u_pos      = np.isfinite(u[i]) and u[i] > 0
        if th_missing and u_pos:
            th[i] = u[i] * 10**r
            n_th += 1

    if verbose and (n_u or n_th):
        print(f'  hpest: estimated {n_u:,} U and {n_th:,} Th from Th/U ratio')


def _predict_rb_from_k(rb, k, verbose):
    """In-place: predict missing Rb from K (wt%) using K/Rb ≈ 300."""
    missing_rb = ~np.isfinite(rb) | (rb <= 0)
    k_pos      = np.isfinite(k) & (k > 0)
    mask       = missing_rb & k_pos
    if mask.any():
        k_ppm = k[mask] * 1e4                   # wt% → ppm
        rb[mask] = k_ppm / _K_RB_RATIO
        if verbose:
            print(f'  hpest: estimated {mask.sum():,} Rb from K (K/Rb = {_K_RB_RATIO:.0f})')


def _predict_sm_from_nd(sm, nd, verbose):
    """In-place: predict missing Sm from Nd using Sm/Nd ≈ 0.22."""
    missing_sm = ~np.isfinite(sm) | (sm <= 0)
    nd_pos     = np.isfinite(nd) & (nd > 0)
    mask       = missing_sm & nd_pos
    if mask.any():
        sm[mask] = nd[mask] * _SM_ND_RATIO
        if verbose:
            print(f'  hpest: estimated {mask.sum():,} Sm from Nd (Sm/Nd = {_SM_ND_RATIO})')


# ---------------------------------------------------------------------------
# U/Th-only estimation for heat production (for backward compatibility)
# ---------------------------------------------------------------------------

def hpest_u_th(data):
    """Estimate missing U or Th from the observed Th/U ratio by rock type.

    Adds ``u_est``, ``th_est``, ``ratio_th_u``, and ``heat_production_est``
    columns.  The main ``hpest()`` function incorporates this logic via
    ``estimate_missing=True``; this function is kept for diagnostic use.

    Parameters
    ----------
    data : pandas.DataFrame
        Requires 'u_ppm', 'th_ppm', 'sio2', 'rock_group', 'rock_origin'.

    Returns
    -------
    pandas.DataFrame
        Modified copy of *data*.
    """
    data = data.copy()
    n    = len(data)

    u_ppm  = _col_or_nan(data, 'u_ppm')
    th_ppm = _col_or_nan(data, 'th_ppm')
    sio2   = _col_or_nan(data, 'sio2')

    has_u  = u_ppm  > 0
    has_th = th_ppm > 0

    data['u_est']      = np.nan
    data['th_est']     = np.nan
    data['ratio_th_u'] = np.nan

    both = has_u & has_th
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = np.where(both, np.log10(th_ppm / u_ppm), np.nan)
    data['ratio_th_u'] = ratio

    grp = (data['rock_group'].str.lower()
           if 'rock_group' in data.columns else pd.Series([''] * n))
    org = (data['rock_origin'].str.lower()
           if 'rock_origin' in data.columns else pd.Series([''] * n))

    ign_mask = ((grp == 'igneous') |
                org.isin(['metaigneous', 'metavolcanic', 'metaplutonic'])).to_numpy()
    sed_mask = ((grp == 'sedimentary') | (org == 'metasedimentary')).to_numpy()

    def _median_ratio(mask):
        r = ratio[mask & both]
        r = r[np.isfinite(r)]
        return float(np.median(r)) if len(r) else np.nan

    med_ign  = _median_ratio(ign_mask)
    print(f'  hpest_u_th: igneous median Th/U = {10**med_ign:.3f}  (log10 = {med_ign:.3f})')

    i_mask = ign_mask & has_th & ~has_u
    data.loc[i_mask, 'u_est'] = th_ppm[i_mask] / 10**med_ign
    print(f'  hpest_u_th: estimated {i_mask.sum():,} U values (igneous)')

    i_mask = ign_mask & has_u & ~has_th
    data.loc[i_mask, 'th_est'] = 10**med_ign * u_ppm[i_mask]
    print(f'  hpest_u_th: estimated {i_mask.sum():,} Th values (igneous)')

    sed_low  = sed_mask & (sio2 <= 30)
    sed_high = sed_mask & (sio2 >= 50)
    sed_mid  = sed_mask & (sio2 > 30) & (sio2 < 50)

    med_low  = _median_ratio(sed_low)
    med_high = _median_ratio(sed_high)

    slope     = (med_high - med_low) / 20.0
    intercept = med_low - 30.0 * slope

    for seg_mask, seg_med, label in [
        (sed_low,  med_low,  'sed SiO2 ≤ 30'),
        (sed_high, med_high, 'sed SiO2 ≥ 50'),
    ]:
        m = seg_mask & has_th & ~has_u
        data.loc[m, 'u_est']  = th_ppm[m] / 10**seg_med
        print(f'  hpest_u_th: estimated {m.sum():,} U values ({label})')

        m = seg_mask & has_u & ~has_th
        data.loc[m, 'th_est'] = 10**seg_med * u_ppm[m]
        print(f'  hpest_u_th: estimated {m.sum():,} Th values ({label})')

    med_mid = slope * sio2 + intercept
    m = sed_mid & has_th & ~has_u
    data.loc[m, 'u_est']  = th_ppm[m] / 10**med_mid[m]
    print(f'  hpest_u_th: estimated {m.sum():,} U values (sed 30–50 SiO2)')

    m = sed_mid & has_u & ~has_th
    data.loc[m, 'th_est'] = 10**med_mid[m] * u_ppm[m]
    print(f'  hpest_u_th: estimated {m.sum():,} Th values (sed 30–50 SiO2)')

    # HP using estimated values
    k2o = get_k2o(data)
    rho = (_col_or_nan(data, 'density_model')
           if 'density_model' in data.columns else np.full(n, np.nan))

    hp_est = np.full(n, np.nan)

    m = (k2o > 0) & has_th & ~has_u
    if m.any():
        hp_est[m] = compute_hp(k2o[m],
                               data.loc[m, 'u_est'].to_numpy(float),
                               th_ppm[m])

    m = (k2o > 0) & has_u & ~has_th
    if m.any():
        hp_est[m] = compute_hp(k2o[m],
                               u_ppm[m],
                               data.loc[m, 'th_est'].to_numpy(float))

    data['heat_production_est'] = hp_est * rho
    return data


# ---------------------------------------------------------------------------
# Basalt liquidus temperature
# ---------------------------------------------------------------------------

def basalt_liquidus(data):
    """Estimate liquidus temperature of basaltic melts from MgO content.

    Applies the empirical relation of Niu et al. (2002)::

        T (°C) = 1026 × exp(0.01894 × MgO)

    Valid for SiO2 < 60 wt%.  Returns NaN for other compositions.

    Parameters
    ----------
    data : pandas.DataFrame
        Requires ``mgo`` (wt%); ``sio2`` applies the SiO2 < 60 filter.

    Returns
    -------
    numpy.ndarray
        Liquidus temperature in °C.

    References
    ----------
    Niu, Y. et al. (2002) Mineral chemistry, whole-rock compositions and
    petrogenesis of ODP Leg 176 gabbros. Proc. Ocean Drill. Prog. Sci.
    Results 176, 1–60.
    """
    mgo  = _col_or_nan(data, 'mgo')
    sio2 = (_col_or_nan(data, 'sio2')
            if 'sio2' in data.columns else np.zeros(len(data)))

    T = 1026.0 * np.exp(0.01894 * mgo)
    T = np.where((sio2 < 60) & np.isfinite(mgo), T, np.nan)
    return T
