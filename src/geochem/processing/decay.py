"""
Radioactive decay correction for elemental concentrations.

``decaycorrect(el, C, age, ref='r88')``
    Given present-day concentration(s) ``C`` of element ``el``, return the
    concentration at time ``age`` Ma in the past by accounting for the
    radioactive decay of that element's radioactive isotopes.

The correction inverts the decay equation:

    C_present = C_stable + Σ_i  abundance_i × C_past × exp(−λ_i × age)

Solving for C_past:

    C_past = C_stable / (1 − Σ abundance_i) is *not* used — instead the
    standard form is:

    C_past = C_stable + Σ_i  abundance_i × C_present × exp(+λ_i × age)

where ``C_stable = (1 − Σ abundance_i) × C_present`` is the non-radioactive
fraction of the element that does not change with time.

Available reference sets (``ref`` parameter):

    ``'hg17'``  Hasterok & Gard (in-house)
    ``'r88'``   Rybach (1988)
    ``'ts14'``  Turcotte & Schubert (2014)
    ``'d12'``   Dye (2012)

Only the five heat-producing parent elements are correctable:
K, Rb, Sm, Th, U.  All other elements are returned unchanged.
"""

from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# Isotope constants  (adapted from radconst.csv)
# Each entry: {'half_life': Ma, 'abundance': fraction of natural element}
# ---------------------------------------------------------------------------

_RADCONST: dict[str, dict[str, list[dict]]] = {
    'hg17': {
        'k':  [{'half_life': 1248.0,     'abundance': 1.170e-4}],
        'rb': [{'half_life': 49700.0,    'abundance': 0.27830}],
        'sm': [{'half_life': 7_000_000., 'abundance': 0.11240}],
        'th': [{'half_life': 14000.0,    'abundance': 1.00000}],
        'u':  [{'half_life': 704.0,      'abundance': 7.204e-3},   # ²³⁵U
               {'half_life': 4468.0,     'abundance': 0.992742}],  # ²³⁸U
    },
    'r88': {
        'k':  [{'half_life': 1230.0,     'abundance': 1.170e-4}],
        'rb': [{'half_life': 49700.0,    'abundance': 1.00000}],
        'sm': [{'half_life': 7_000_000., 'abundance': 1.00000}],
        'th': [{'half_life': 13900.0,    'abundance': 1.00000}],
        'u':  [{'half_life': 713.0,      'abundance': 7.110e-3},   # ²³⁵U
               {'half_life': 4510.0,     'abundance': 0.99280}],   # ²³⁸U
    },
    'ts14': {
        'k':  [{'half_life': 1250.0,     'abundance': 1.190e-4}],
        'rb': [{'half_life': 49700.0,    'abundance': 1.00000}],
        'sm': [{'half_life': 7_000_000., 'abundance': 1.00000}],
        'th': [{'half_life': 14000.0,    'abundance': 1.00000}],
        'u':  [{'half_life': 704.0,      'abundance': 7.100e-3},   # ²³⁵U
               {'half_life': 4470.0,     'abundance': 0.99280}],   # ²³⁸U
    },
    'd12': {
        'k':  [{'half_life': 1280.0,     'abundance': 1.170e-4}],
        'rb': [{'half_life': 49700.0,    'abundance': 1.00000}],
        'sm': [{'half_life': 7_000_000., 'abundance': 1.00000}],
        'th': [{'half_life': 14000.0,    'abundance': 1.00000}],
        'u':  [{'half_life': 704.0,      'abundance': 7.204e-3},   # ²³⁵U
               {'half_life': 4460.0,     'abundance': 0.992796}],  # ²³⁸U
    },
}

# Elements correctable by this module (case-insensitive, _ppm suffix stripped)
DECAY_ELEMENTS: frozenset[str] = frozenset({'k', 'rb', 'sm', 'th', 'u'})


def decaycorrect(el: str, C, age, ref: str = 'r88') -> np.ndarray:
    """Correct element concentration(s) to their values at *age* Ma in the past.

    Parameters
    ----------
    el : str
        Element name.  Case-insensitive; a trailing ``'_ppm'`` suffix is
        stripped (e.g. ``'th_ppm'`` is treated the same as ``'Th'``).
        Correctable elements: K, Rb, Sm, Th, U.  Any other element is
        returned unchanged.
    C : array-like
        Present-day concentration(s) in ppm (or any consistent unit).
        Must be non-negative; negative / NaN values are passed through
        unchanged.
    age : float or array-like
        Time in the past in Ma.  Scalar is broadcast to the length of *C*.
    ref : {'r88', 'hg17', 'ts14', 'd12'}
        Reference for isotope constants (default ``'r88'``, Rybach 1988).

    Returns
    -------
    numpy.ndarray
        Corrected concentrations in the same units as *C*.

    Examples
    --------
    >>> decaycorrect('Th', 10.0, 2500)   # Th at 2500 Ma
    >>> decaycorrect('u_ppm', data['u_ppm'].to_numpy(), data['avg_age'].to_numpy())

    Notes
    -----
    The correction:

        C_past = C_stable + Σᵢ abundance_i × C × exp(λᵢ × age)

    where ``C_stable = (1 − Σ abundance_i) × C`` and ``λᵢ = ln2 / t½_i``.

    For Rb and Sm with ``ref='r88'``, ``'ts14'``, or ``'d12'``, ``abundance=1``
    so ``C_stable = 0`` — all atoms are assumed to have been radioactive
    (a conservative approximation for these minor contributors).

    References
    ----------
    Rybach (1988) *Determination of heat production rate*, in Haenel et al.,
    Handbook of Terrestrial Heat-Flow Density Determination.
    """
    ref_key = ref.strip().lower()
    if ref_key not in _RADCONST:
        raise ValueError(
            f"Unknown ref '{ref}'.  Choose from: {sorted(_RADCONST)}."
        )

    # Strip _ppm suffix, normalise to lowercase
    el_key = el.strip().lower()
    if el_key.endswith('_ppm'):
        el_key = el_key[:-4]

    C = np.asarray(C, dtype=float).ravel()

    if el_key not in DECAY_ELEMENTS:
        return C.copy()

    isotopes = _RADCONST[ref_key][el_key]

    age = np.asarray(age, dtype=float).ravel()
    if age.size == 1:
        age = np.broadcast_to(age, C.shape)

    # Stable fraction of the element (does not change with time)
    total_abundance = sum(iso['abundance'] for iso in isotopes)
    C_stable = (1.0 - total_abundance) * C

    # Sum radiogenic contributions decayed back to time = age
    C_past = C_stable.copy()
    for iso in isotopes:
        lam = np.log(2.0) / iso['half_life']  # 1/Ma
        C_past = C_past + iso['abundance'] * C * np.exp(lam * age)

    # Preserve NaN / negative input values unchanged
    bad = ~np.isfinite(C) | (C < 0)
    C_past[bad] = C[bad]

    return C_past
