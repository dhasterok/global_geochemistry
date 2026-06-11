"""
Geochemical classification indices for igneous and sedimentary rocks.

All input oxide columns are assumed to be in wt% with iron as FeO_tot.
Below-detection-limit values (≤ 0) are treated as zero.

References:
    Frost et al. (2001) A geochemical classification for granitic rocks.
    Journal of Petrology 42, 2033–2048.

    Debon & LeFort (1983) A chemical-mineralogical classification of common
    plutonic rocks. Trans. R. Soc. Edinb. Earth Sci. 73, 135–149.

    De La Roche et al. (1980) A classification of volcanic and plutonic rocks
    using R1–R2 diagram. Earth and Planetary Science Letters 48, 69–79.

    Price & Velbel (2003) Chemical weathering indices applied to weathering
    profiles developed on heterogeneous felsic metamorphic parent rocks.
    Chemical Geology 202, 397–416.

    O'Neill, H.St.C. (2016) The smoothness and shapes of chondrite-normalised
    rare earth element patterns in basalts. Journal of Petrology 57, 1463–1508.
"""

import numpy as np
import pandas as pd
from src.utils.molecular import MolecularWeightCalculator

_mwc = MolecularWeightCalculator()

# Precomputed molecular weights for the common major oxides
_MW = {ox: _mwc.molecular_weight(ox) for ox in
       ['SiO2', 'TiO2', 'Al2O3', 'FeO', 'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5']}


def geochem_index(data, inplace=False):
    """Compute standard geochemical classification indices.

    Requires normalized oxide columns (lowercase): sio2, tio2, al2o3,
    feo_tot, mgo, cao, na2o, k2o, p2o5.  Missing columns are treated as zero.

    Parameters
    ----------
    data : pandas.DataFrame
    inplace : bool

    Returns
    -------
    pandas.DataFrame
        Input data with the following columns added:

        Mg_number   Mg/(Mg + 0.85×Fe) molar — Frost et al. (2001)
        Fe_number   FeO/(FeO + MgO) by weight — Frost et al. (2001)
        MALI        Modified Alkali-Lime Index: Na2O + K2O − CaO — Frost et al.
        CAI         Calcic-Alkalic Index (deviation from calc-alkalic line)
        ASI         Aluminium Saturation Index — Frost et al. (2001)
        AI          Agpaitic Index: Al2O3 / (Na2O + K2O) molar
        maficity    Mg + Fe + Ti molar — Debon & LeFort (1983)
        spar        Feldspar index — Debon & LeFort (1983)
        qtzindex    Quartz index — Debon & LeFort (1983)
        CIA         Chemical Index of Alteration — Price & Velbel (2003)
        WIP         Weathering Index of Parker — Price & Velbel (2003)
        CPA         Chemical Proxy of Alteration — Buggle et al. (2011)
        R1, R2      Roche et al. (1980) discriminants
    """
    if not inplace:
        data = data.copy()

    def _get(col):
        if col in data.columns:
            v = data[col].copy()
        else:
            v = pd.Series(0.0, index=data.index)
        v = v.where(v > 0, 0.0)    # BDL → 0
        return v

    sio2  = _get('sio2')
    tio2  = _get('tio2')
    al2o3 = _get('al2o3')
    feo   = _get('feo_tot')
    mgo   = _get('mgo')
    cao   = _get('cao')
    na2o  = _get('na2o')
    k2o   = _get('k2o')
    p2o5  = _get('p2o5')

    # Mole fractions (wt% / molecular weight)
    n_sio2  = sio2  / _MW['SiO2']
    n_tio2  = tio2  / _MW['TiO2']
    n_al2o3 = al2o3 / _MW['Al2O3']
    n_feo   = feo   / _MW['FeO']
    n_mgo   = mgo   / _MW['MgO']
    n_cao   = cao   / _MW['CaO']
    n_na2o  = na2o  / _MW['Na2O']
    n_k2o   = k2o   / _MW['K2O']
    n_p2o5  = p2o5  / _MW['P2O5']

    # CaO adjusted for apatite contribution (remove 10/3 P2O5 moles of CaO)
    n_ca_noap = (n_cao - (10.0 / 3.0) * n_p2o5).clip(lower=0)

    def _safe(series):
        """Replace inf with NaN."""
        return series.replace([np.inf, -np.inf], np.nan)

    # Mg number: Mg/(Mg + 0.85 Fe) molar
    data['Mg_number'] = _safe(n_mgo / (n_mgo + 0.85 * n_feo))

    # Fe number: FeO/(FeO + MgO) by weight — Frost et al. (2001)
    data['Fe_number'] = _safe(feo / (feo + mgo))

    # MALI: Modified Alkali-Lime Index
    data['MALI'] = na2o + k2o - cao

    # CAI: Calcic-Alkalic Index (deviation from dividing line)
    data['CAI'] = data['MALI'] - (-44.72 + 1.094 * sio2 - 0.00527 * sio2 ** 2)

    # ASI: Aluminium Saturation Index — Frost et al. (2001)
    data['ASI'] = _safe(n_al2o3 / (n_ca_noap + n_na2o + n_k2o))

    # AI: Agpaitic Index
    data['AI'] = _safe(n_al2o3 / (n_na2o + n_k2o))

    # Maficity — Debon & LeFort (1983)
    data['maficity'] = n_mgo + n_feo + n_tio2

    # Feldspar index — Debon & LeFort (1983)
    data['spar'] = 2 * n_k2o - (2 * n_na2o + n_cao - (10.0 / 3.0) * n_p2o5)

    # Quartz index — Debon & LeFort (1983)
    data['qtzindex'] = n_sio2 / 3.0 - (2 * n_na2o + 2 * n_k2o + n_ca_noap / 3.0)

    # CIA — Price & Velbel (2003)
    data['CIA'] = _safe(
        100 * n_al2o3 / (n_al2o3 + n_ca_noap + n_na2o + n_k2o)
    )

    # WIP — Price & Velbel (2003)
    data['WIP'] = 100 * (2 * n_na2o / 0.35 + n_mgo / 0.9 + 8 * n_k2o + n_cao / 0.7)

    # CPA — Buggle et al. (2011)
    data['CPA'] = _safe(100 * n_al2o3 / (n_al2o3 + n_na2o))

    # R1 & R2 — De La Roche et al. (1980)
    data['R1'] = 4 * n_sio2 - 22 * (n_na2o + n_k2o) - 2 * (n_feo + n_tio2)
    data['R2'] = 6 * n_cao + 2 * n_mgo + 2 * n_al2o3

    return data


# ---------------------------------------------------------------------------
# REE lambda coefficients (O'Neill 2016)
# ---------------------------------------------------------------------------

# REE column names expected in the DataFrame
_REE_COLS = [
    'la_ppm', 'ce_ppm', 'pr_ppm', 'nd_ppm', 'sm_ppm', 'eu_ppm', 'gd_ppm',
    'tb_ppm', 'dy_ppm', 'ho_ppm', 'er_ppm', 'tm_ppm', 'yb_ppm', 'lu_ppm',
]
_EU_IDX = 5  # index of Eu in the 14-element list

# REE ionic radii in Å — Shannon (1976), 8-coordinate
_REE_R = np.array([
    1.160, 1.143, 1.126, 1.109, 1.079, 1.066, 1.053,
    1.040, 1.027, 1.015, 1.004, 0.994, 0.985, 0.977,
])

# CI chondrite normalisation — O'Neill (2016) Table 1
_CI = np.array([
    0.2472, 0.6308, 0.0950, 0.4793, 0.15419, 0.0592, 0.2059,
    0.0375, 0.2540, 0.0554, 0.1645, 0.0258,  0.1684, 0.0251,
])

# Orthonormal polynomial root constants — O'Neill (2016)
_B1 = _REE_R - 1.054769
_B2 = (_REE_R - 1.005327429) * (_REE_R - 1.128236038)
_B3 = (_REE_R - 1.060548105) * (_REE_R - 1.145519887) * (_REE_R - 0.991412204)
_B4 = ((_REE_R - 1.104414973) * (_REE_R - 1.153426708)
       * (_REE_R - 0.984820219) * (_REE_R - 1.030518142))

# Basis matrix (5 × 14): basis[k, i] = k-th basis function at REE i
_BASIS = np.vstack([np.ones(14), _B1, _B2, _B3, _B4])

# Operator sub-matrices (5 × 14): A_all[k, j, i] = r_i^j * Bk_i
_A_ALL = np.stack([
    _REE_R[np.newaxis, :] ** np.arange(5)[:, np.newaxis],                  # A1: r^j
    (_REE_R[np.newaxis, :] ** np.arange(5)[:, np.newaxis]) * _B1,          # A2: r^j * B1
    (_REE_R[np.newaxis, :] ** np.arange(5)[:, np.newaxis]) * _B2,          # A3: r^j * B2
    (_REE_R[np.newaxis, :] ** np.arange(5)[:, np.newaxis]) * _B3,          # A4: r^j * B3
    (_REE_R[np.newaxis, :] ** np.arange(5)[:, np.newaxis]) * _B4,          # A5: r^j * B4
])  # shape (5, 5, 14): [basis_idx, poly_order, ree_idx]

_EU_MASK = np.zeros(14, dtype=bool)
_EU_MASK[_EU_IDX] = True


def lambdaree(data, normalize=True, min_ree=7):
    """Orthonormal polynomial fit to chondrite-normalised REE patterns.

    Computes lambda0–lambda4 coefficients and a chi-squared misfit (MSWD)
    for each sample, following the method of O'Neill (2016).  Eu is excluded
    from the fit because of its susceptibility to redox fractionation.

    Parameters
    ----------
    data : pandas.DataFrame
        Must contain REE columns in ppm (la_ppm … lu_ppm).  Missing or
        non-positive values are treated as NaN.
    normalize : bool
        Divide by CI chondrite before taking log (default True).
    min_ree : int
        Minimum number of non-NaN, non-Eu REE required to attempt a fit
        (default 7).

    Returns
    -------
    pandas.DataFrame
        Input DataFrame with seven new columns:
        lambda0–lambda4, nree (number used), lambda_mswd.

    References
    ----------
    O'Neill, H.St.C. (2016) The smoothness and shapes of chondrite-normalised
    rare earth element patterns in basalts. J. Petrology 57, 1463–1508.
    Original algorithm by Hugh O'Neill (Excel, Dec 2013); MATLAB port by
    D. Hasterok (2020).
    """
    n = len(data)

    # Build (n × 14) sample matrix; set ≤0 values to NaN
    sample = np.full((n, 14), np.nan)
    for j, col in enumerate(_REE_COLS):
        if col in data.columns:
            vals = data[col].to_numpy(float)
            sample[:, j] = np.where(vals <= 0, np.nan, vals)

    with np.errstate(divide='ignore', invalid='ignore'):
        if normalize:
            sample = np.log(sample / _CI[np.newaxis, :])
        else:
            sample = np.log(sample)

    lam   = np.full((n, 5), np.nan)
    nree  = np.full(n, np.nan)
    chi2  = np.full(n, np.nan)

    for k in range(n):
        X   = sample[k]
        ind = ~np.isnan(X) & ~_EU_MASK
        if ind.sum() < min_ree:
            continue

        nree[k] = ind.sum()

        # Operator matrix A (5 × 5): col m = sum over selected REE of _A_ALL[m, :, i]
        # Index in two steps to avoid numpy mixed advanced/basic indexing reordering
        A_mat = np.column_stack(
            [_A_ALL[m][:, ind].sum(axis=1) for m in range(5)]
        )  # (5, 5)

        # Data vector d (5,): d[j] = Σ_i X[i] * r[i]^j  (using A1 = _A_ALL[0])
        d_vec = (_A_ALL[0][:, ind] * X[ind]).sum(axis=1)   # (5,)

        try:
            Q = np.linalg.solve(A_mat, d_vec)
        except np.linalg.LinAlgError:
            continue

        lam[k] = Q

        Xm       = Q @ _BASIS                       # (14,)
        resid    = X[ind] - Xm[ind]
        chi2[k]  = (resid**2).sum() / (1e-4 * (ind.sum() - 4))

    out = data.copy()
    for i, name in enumerate(['lambda0', 'lambda1', 'lambda2', 'lambda3', 'lambda4']):
        out[name] = lam[:, i]
    out['nree']         = nree
    out['lambda_mswd']  = chi2
    return out
