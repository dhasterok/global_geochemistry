"""
CIPW normative mineralogy.

Implements the sequential mineral assignment of Cross, Iddings, Pirsson &
Washington.  All mole calculations follow Verma et al. (2002, 2003).  The
iron speciation formula is from Le Maitre (1976).

References:
    Verma, S.P., Torres-Alvarado, I.S., Sotelo-Rodriguez, Z.T. (2002)
    SINCLAS: standard igneous norm and volcanic rock classification system.
    Computers & Geosciences 28, 711–724.

    Le Maitre, R.W. (1976) Some problems of the projection of chemical data
    into mineralogical classifications.  Contrib. Mineral. Petrol. 56, 181–189.

Note on MATLAB bug fixed here:
    cipwnorm.m line 192 used m.Tn (zero at that point) instead of m.Tnp when
    reducing CaO during the Sphene/Rutile step.  This is corrected below.
"""

import numpy as np
import pandas as pd

from src.utils.molecular import MolecularWeightCalculator

_mwc = MolecularWeightCalculator()

# Frequently used molecular weights (g/mol)
_MW = {
    ox: _mwc.molecular_weight(ox)
    for ox in ['SiO2', 'TiO2', 'Al2O3', 'Fe2O3', 'FeO', 'MnO',
               'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'CO2']
}
_Rfe = 2 * _MW['FeO'] / _MW['Fe2O3']   # ≈ 0.8998

# Constitutional (formula) weights of normative minerals (g/mol)
# Source: Verma et al. (2002) Table 2 / minweight.csv
_W = {
    'Qz':     60.0843,      # SiO2
    'C':     101.961276,    # Al2O3
    'Or':    556.663076,    # K2O.Al2O3.6SiO2
    'Ab':    524.446016,    # Na2O.Al2O3.6SiO2
    'An':    278.207276,    # CaO.Al2O3.2SiO2
    'Ne':    284.108816,    # Na2O.Al2O3.2SiO2
    'Lc':    436.494476,    # K2O.Al2O3.4SiO2
    'Kp':    316.325876,    # K2O.Al2O3.2SiO2
    'Ac':    462.004339,    # Na2O.Fe2O3.4SiO2
    'Ns':    122.063239,    # Na2O.SiO2
    'Ks':    154.2803,      # K2O.SiO2
    'Di_Mg': 216.5504,      # CaO.MgO.2SiO2
    'Di_Fe_corr': 176.246,  # CaO + 2SiO2 part (FeO/MnO added at runtime)
    'Wo':    116.1617,      # CaO.SiO2
    'Hy_Mg': 100.3886,      # MgO.SiO2
    'Hy_Fe_corr':  60.0843, # SiO2 part (FeO/MnO added at runtime)
    'Ol_Mg': 140.6931,      # 2MgO.SiO2
    'Ol_Fe_corr':  60.0843, # SiO2 part (2x FeO/MnO added at runtime)
    'Cs':    172.2391,      # 2CaO.SiO2
    'Tn':    196.0275,      # CaO.TiO2.SiO2
    'Pf':    135.9432,      # CaO.TiO2
    'Ap':    328.869191,    # 3CaO.P2O5.(1/3)CaO
    'Nc':    105.988339,    # Na2O.CO2
    'Cc':    100.0868,      # CaO.CO2
    'He':    159.6882,      # Fe2O3
    'Ru':     79.8658,      # TiO2
    'Mt_corr': 159.6882,    # Fe2O3 part (FeO/MnO added at runtime)
    'Il_corr':  79.8658,    # TiO2 part (FeO/MnO added at runtime)
}


def _fe_ratio_lemaitre(sio2, na2o, k2o, is_volcanic):
    """Fe2O3/FeO ratio from Le Maitre (1976) empirical formula."""
    ratio = np.where(
        is_volcanic,
        0.93 - 0.0042 * sio2 - 0.022 * (na2o + k2o),
        0.88 - 0.0016 * sio2 - 0.027 * (na2o + k2o),
    )
    return np.clip(ratio, 0.0, None)


def _prepare(data, is_volcanic, use_measured_fe=True):
    """Anhydrous normalisation and Fe2O3/FeO split.

    Returns arrays (all per 100 g anhydrous):
        n_sio2, n_tio2, n_al2o3, n_feo, n_fe2o3, n_mno,
        n_mgo, n_cao, n_na2o, n_k2o, n_p2o5, n_co2,
        xfeo, xmno, xfer, xmgr
    """
    n = len(data)

    def col(name, default=0.0):
        arr = data[name].to_numpy(float) if name in data.columns else np.full(n, default)
        arr = np.nan_to_num(arr, nan=0.0)
        return np.clip(arr, 0.0, None)   # negative values → 0

    sio2  = col('sio2')
    tio2  = col('tio2')
    al2o3 = col('al2o3')
    mno   = col('mno')
    mgo   = col('mgo')
    cao   = col('cao')
    na2o  = col('na2o')
    k2o   = col('k2o')
    p2o5  = col('p2o5')
    co2   = col('co2')

    # Determine feo and fe2o3 from whatever iron columns are present
    feo_tot = col('feo_tot') if 'feo_tot' in data.columns else col('feo')
    feo_raw = col('feo')
    fe2o3_raw = col('fe2o3')

    has_both = (feo_raw > 0) & (fe2o3_raw > 0)
    feot = np.where(has_both, feo_raw + _Rfe * fe2o3_raw, feo_tot)

    # Fe ratio: use measured when both present, else Le Maitre (1976)
    fe_ratio = _fe_ratio_lemaitre(sio2, na2o, k2o, is_volcanic)
    if use_measured_fe:
        with np.errstate(invalid='ignore', divide='ignore'):
            raw_ratio = np.where(feo_raw > 0, fe2o3_raw / feo_raw, 0.0)
        fe_ratio = np.where(has_both, raw_ratio, fe_ratio)

    feo   = feot / (1.0 + _Rfe * fe_ratio)
    fe2o3 = (feot - feo) / _Rfe

    # Anhydrous normalisation (oxides 1–11, excluding LOI and CO2)
    oxsum = sio2 + tio2 + al2o3 + feo + fe2o3 + mno + mgo + cao + na2o + k2o + p2o5
    oxsum = np.where(oxsum <= 0, np.nan, oxsum)
    scale = 100.0 / oxsum

    def norm(x): return x * scale

    sio2  = norm(sio2);   tio2 = norm(tio2);  al2o3 = norm(al2o3)
    feo   = norm(feo);  fe2o3 = norm(fe2o3);  mno   = norm(mno)
    mgo   = norm(mgo);  cao   = norm(cao);    na2o  = norm(na2o)
    k2o   = norm(k2o);  p2o5  = norm(p2o5)

    # Mole fractions (per 100 g)
    n_sio2  = sio2  / _MW['SiO2']
    n_tio2  = tio2  / _MW['TiO2']
    n_al2o3 = al2o3 / _MW['Al2O3']
    n_fe2o3 = fe2o3 / _MW['Fe2O3']
    n_mno   = mno   / _MW['MnO']
    n_mgo   = mgo   / _MW['MgO']
    n_cao   = cao   / _MW['CaO']
    n_na2o  = na2o  / _MW['Na2O']
    n_k2o   = k2o   / _MW['K2O']
    n_p2o5  = p2o5  / _MW['P2O5']
    n_co2   = co2   / _MW['CO2']   # not normalised; used directly

    # n_feo lumps FeO + MnO (Verma convention)
    n_feo = feo / _MW['FeO'] + n_mno

    # FeO/MnO partition within the combined Fe+Mn mole pool.
    # These fractions are fixed at the start (MnO never oxidised separately).
    xfeo = np.where(n_feo > 0, feo / _MW['FeO'] / n_feo, 0.0)
    xmno = np.where(n_feo > 0, n_mno / n_feo, 0.0)

    # xfer / xmgr (Mg vs Fe+Mn partition) are computed in cipwnorm AFTER
    # ilmenite and magnetite steps consume some of the n_feo pool.

    return (n_sio2, n_tio2, n_al2o3, n_feo, n_fe2o3, n_mno,
            n_mgo, n_cao, n_na2o, n_k2o, n_p2o5, n_co2,
            xfeo, xmno)


def cipwnorm(data, is_volcanic=None, co2_mode='ignore'):
    """Compute CIPW normative mineralogy.

    Parameters
    ----------
    data : pandas.DataFrame
        Whole-rock major oxide data (wt%).  Must contain at minimum 'sio2',
        'feo_tot' (or 'feo'/'fe2o3'), 'mgo', 'cao', 'na2o', 'k2o'.
    is_volcanic : array-like of bool, optional
        True for volcanic rocks, False for plutonic.  Used for Le Maitre
        (1976) Fe-ratio calibration.  Defaults to True for all samples.
    co2_mode : {'ignore', 'cancrinite', 'calcite'}
        How to handle CO2.  'ignore' (default) sets CO2 to 0.  'cancrinite'
        assigns CO2 to Na2O first, then CaO.  'calcite' assigns CO2 to CaO.

    Returns
    -------
    pandas.DataFrame
        Normative mineral proportions (wt%), one column per mineral, plus
        derived indices.  Same index as *data*.

    Notes
    -----
    Mineral abbreviations follow Kretz (1983) / Siivola & Schmid (2007).
    """
    n = len(data)
    if is_volcanic is None:
        is_volcanic = np.ones(n, dtype=bool)
    else:
        is_volcanic = np.asarray(is_volcanic, dtype=bool)

    (n_sio2, n_tio2, n_al2o3, n_feo, n_fe2o3, n_mno,
     n_mgo, n_cao, n_na2o, n_k2o, n_p2o5, n_co2,
     xfeo, xmno) = _prepare(data, is_volcanic)

    # Working copies (moles, consumed as minerals are assigned)
    q_sio2  = n_sio2.copy()
    q_tio2  = n_tio2.copy()
    q_al2o3 = n_al2o3.copy()
    q_feo   = n_feo.copy()
    q_fe2o3 = n_fe2o3.copy()
    q_cao   = n_cao.copy()
    q_na2o  = n_na2o.copy()
    q_k2o   = n_k2o.copy()
    q_p2o5  = n_p2o5.copy()
    q_co2   = n_co2.copy() if co2_mode != 'ignore' else np.zeros(n)

    # Provisional mineral moles (intermediate values)
    # Final minerals are named without suffix; provisionals have _p suffix
    def _alloc():
        return np.zeros(n)

    Ap_p = _alloc()    # apatite provisional
    Nc   = _alloc()    # cancrinite (sodium carbonate)
    Cc   = _alloc()    # calcite
    Il   = _alloc()    # ilmenite
    Orp  = _alloc()    # orthoclase provisional
    Ks   = _alloc()    # potassium metasilicate
    Abp  = _alloc()    # albite provisional
    Ac   = _alloc()    # acmite
    Ns   = _alloc()    # sodium metasilicate
    An   = _alloc()    # anorthite
    Cor  = _alloc()    # corundum (excess Al2O3)
    Tnp  = _alloc()    # sphene provisional
    Ru   = _alloc()    # rutile
    Mt   = _alloc()    # magnetite
    Hm   = _alloc()    # hematite
    Dip  = _alloc()    # diopside provisional
    Wop  = _alloc()    # wollastonite provisional
    Hyp  = _alloc()    # hypersthene provisional

    # ------------------------------------------------------------------ #
    # Step 1: Apatite  (uses P2O5 and CaO)
    # Ca5(PO4)3(OH,F,Cl) → 3CaO.P2O5.(1/3)CaO  → ratio CaO:P2O5 = 10/3
    # ------------------------------------------------------------------ #
    ind = q_cao >= (10/3) * q_p2o5
    Ap_p[ ind] = q_p2o5[ind]
    q_cao[ ind] -= (10/3) * Ap_p[ind]
    q_p2o5[ind]  = 0.0

    Ap_p[~ind] = (3/10) * q_cao[~ind]
    q_p2o5[~ind] -= Ap_p[~ind]
    q_cao[~ind]  = 0.0

    free_p2o5 = q_p2o5.copy()    # any P2O5 not in apatite

    # ------------------------------------------------------------------ #
    # Step 2: CO2 minerals (cancrinite or calcite; default: skip)
    # ------------------------------------------------------------------ #
    if co2_mode == 'cancrinite':
        ind = q_na2o >= q_co2
        Nc[ ind] = q_co2[ind]
        q_na2o[ ind] -= Nc[ind]
        q_co2[  ind]  = 0.0

        Nc[~ind] = q_na2o[~ind]
        q_co2[~ind] -= Nc[~ind]
        q_na2o[~ind]  = 0.0

    elif co2_mode == 'calcite':
        ind = q_cao >= q_co2
        Cc[ ind] = q_co2[ind]
        q_cao[ ind] -= Cc[ind]
        q_co2[ ind]  = 0.0

        Cc[~ind] = q_cao[~ind]
        q_co2[~ind] -= Cc[~ind]
        q_cao[~ind]  = 0.0

    free_co2 = q_co2.copy()

    # ------------------------------------------------------------------ #
    # Step 3: Ilmenite  (uses TiO2 and FeO+MnO)
    # ------------------------------------------------------------------ #
    ind = q_feo >= q_tio2
    Il[ ind] = q_tio2[ind]
    q_feo[ ind] -= Il[ind]
    q_tio2[ind]  = 0.0

    Il[~ind] = q_feo[~ind]
    q_tio2[~ind] -= Il[~ind]
    q_feo[~ind]   = 0.0

    # ------------------------------------------------------------------ #
    # Step 4: Orthoclase provisional  (uses K2O and Al2O3)
    # ------------------------------------------------------------------ #
    ind = q_al2o3 >= q_k2o
    Orp[ ind] = q_k2o[ind]
    q_al2o3[ind] -= Orp[ind]
    q_k2o[  ind]  = 0.0

    Orp[~ind] = q_al2o3[~ind]
    q_k2o[~ind] -= Orp[~ind]
    q_al2o3[~ind] = 0.0

    Ks = q_k2o.copy()     # excess K2O → potassium metasilicate

    Y = 6 * Orp + Ks      # running SiO2 tally (moles needed by minerals so far)

    # ------------------------------------------------------------------ #
    # Step 5: Albite provisional  (uses Na2O and Al2O3)
    # ------------------------------------------------------------------ #
    ind = q_al2o3 >= q_na2o
    Abp[ ind] = q_na2o[ind]
    q_al2o3[ind] -= Abp[ind]
    q_na2o[ ind]  = 0.0

    Abp[~ind] = q_al2o3[~ind]
    q_na2o[~ind] -= Abp[~ind]
    q_al2o3[~ind] = 0.0

    Y += 6 * Abp

    # ------------------------------------------------------------------ #
    # Step 6: Acmite / sodium metasilicate  (uses excess Na2O and Fe2O3)
    # ------------------------------------------------------------------ #
    ind = q_na2o >= q_fe2o3
    Ac[ ind] = q_fe2o3[ind]
    q_na2o[ ind] -= Ac[ind]
    q_fe2o3[ind]  = 0.0

    Ac[~ind] = q_na2o[~ind]
    q_fe2o3[~ind] -= Ac[~ind]
    q_na2o[~ind]  = 0.0

    Ns = q_na2o.copy()   # residual Na2O → sodium metasilicate
    Y += 4 * Ac + Ns

    # ------------------------------------------------------------------ #
    # Step 7: Anorthite / corundum  (uses CaO and Al2O3)
    # ------------------------------------------------------------------ #
    ind = q_al2o3 >= q_cao
    An[ ind] = q_cao[ind]
    q_al2o3[ind] -= An[ind]
    q_cao[  ind]  = 0.0

    An[~ind] = q_al2o3[~ind]
    q_cao[~ind] -= An[~ind]
    q_al2o3[~ind] = 0.0

    Cor = q_al2o3.copy()   # residual Al2O3 → corundum
    Y += 2 * An

    # ------------------------------------------------------------------ #
    # Step 8: Sphene provisional / rutile  (uses TiO2 and CaO)
    # Bug fix: MATLAB used m.Tn (=0) instead of m.Tnp here.
    # ------------------------------------------------------------------ #
    ind = q_cao >= q_tio2
    Tnp[ ind] = q_tio2[ind]
    q_cao[ ind] -= Tnp[ind]    # fixed: was -= m.Tn (always 0 in MATLAB)
    q_tio2[ind]  = 0.0

    Tnp[~ind] = q_cao[~ind]
    q_tio2[~ind] -= Tnp[~ind]
    q_cao[~ind]  = 0.0

    Ru = q_tio2.copy()   # residual TiO2 → rutile
    Y += Tnp

    # ------------------------------------------------------------------ #
    # Step 9: Magnetite / hematite  (uses Fe2O3 and FeO+MnO)
    # ------------------------------------------------------------------ #
    ind = q_fe2o3 >= q_feo
    Mt[ ind] = q_feo[ind]
    q_fe2o3[ind] -= Mt[ind]
    q_feo[  ind]  = 0.0

    Mt[~ind] = q_fe2o3[~ind]
    q_feo[~ind] -= Mt[~ind]
    q_fe2o3[~ind] = 0.0

    Hm = q_fe2o3.copy()   # residual Fe2O3 → hematite

    # Mg vs Fe+Mn partition — computed here using REMAINING q_feo (after Il + Mt),
    # matching MATLAB cipwnorm.m lines 217–220.
    n_femg_remaining = q_feo + n_mgo
    xfer = np.where(n_femg_remaining > 0, q_feo / n_femg_remaining, 0.0)
    xmgr = np.where(n_femg_remaining > 0, n_mgo / n_femg_remaining, 0.0)

    # ------------------------------------------------------------------ #
    # Step 10: Diopside / wollastonite / hypersthene
    # ------------------------------------------------------------------ #
    n_femg = n_femg_remaining

    ind = q_cao >= n_femg
    Dip[ ind] = n_femg[ind]
    q_cao[ ind] -= Dip[ind]
    n_femg[ind]  = 0.0

    Dip[~ind] = q_cao[~ind]
    n_femg[~ind] -= Dip[~ind]
    q_cao[~ind]  = 0.0

    Wop = q_cao.copy()      # residual CaO → wollastonite
    Hyp = n_femg.copy()     # residual (Fe,Mg)O → hypersthene

    Y += 2 * Dip + Wop + Hyp

    # ------------------------------------------------------------------ #
    # Step 11: Quartz / undersaturation deficit
    # ------------------------------------------------------------------ #
    Q = np.maximum(q_sio2 - Y, 0.0)
    D = np.maximum(Y - q_sio2, 0.0)
    unsat = D > 0.0

    # ------------------------------------------------------------------ #
    # Step 12: Olivine / hypersthene adjustment
    # ------------------------------------------------------------------ #
    Olp  = _alloc()
    Hy   = _alloc()
    D1   = _alloc()

    ind = D < Hyp / 2
    Olp[ unsat &  ind] = D[unsat & ind]
    Hy[  unsat &  ind] = (Hyp - 2*D)[unsat & ind]
    D1[  unsat &  ind] = 0.0

    Olp[ unsat & ~ind] = Hyp[unsat & ~ind] / 2
    Hy[  unsat & ~ind] = 0.0
    D1[  unsat & ~ind] = (D - Hyp/2)[unsat & ~ind]

    # Saturated: no olivine, keep full Hy
    Hy[~unsat] = Hyp[~unsat]

    unsat = D1 > 0.0

    # ------------------------------------------------------------------ #
    # Step 13: Sphene / perovskite
    # ------------------------------------------------------------------ #
    Pf = _alloc()
    Tn = _alloc()
    D2 = _alloc()

    ind = D1 < Tnp
    Pf[ unsat &  ind] = D1[unsat & ind]
    Tn[ unsat &  ind] = (Tnp - D1)[unsat & ind]
    D2[ unsat &  ind] = 0.0

    Pf[ unsat & ~ind] = Tnp[unsat & ~ind]
    Tn[ unsat & ~ind] = 0.0
    D2[ unsat & ~ind] = (D1 - Tnp)[unsat & ~ind]

    unsat = D2 > 0.0

    # ------------------------------------------------------------------ #
    # Step 14: Nepheline / albite
    # ------------------------------------------------------------------ #
    Ne   = _alloc()
    Abpp = _alloc()
    Ab   = _alloc()
    D3   = _alloc()

    ind = D2 < 4 * Abp
    Ne[   unsat &  ind] = D2[unsat & ind] / 4
    Abpp[ unsat &  ind] = (Abp - D2/4)[unsat & ind]
    D3[   unsat &  ind] = 0.0

    Ne[   unsat & ~ind] = Abp[unsat & ~ind]
    Abpp[ unsat & ~ind] = 0.0
    D3[   unsat & ~ind] = (D2 - 4*Abp)[unsat & ~ind]

    Ab[ unsat] = Abpp[unsat]
    Ab[~unsat] = Abp[~unsat]

    unsat = D3 > 0.0

    # ------------------------------------------------------------------ #
    # Step 15: Leucite / orthoclase
    # ------------------------------------------------------------------ #
    Lcp  = _alloc()
    Orpp = _alloc()
    Or   = _alloc()
    D4   = _alloc()

    ind = D3 < 2 * Orp
    Lcp[  unsat &  ind] = D3[unsat & ind] / 2
    Orpp[ unsat &  ind] = (Orp - D3/2)[unsat & ind]
    D4[   unsat &  ind] = 0.0

    Lcp[  unsat & ~ind] = Orp[unsat & ~ind]
    Orpp[ unsat & ~ind] = 0.0
    D4[   unsat & ~ind] = (D3 - 2*Orp)[unsat & ~ind]

    Or[ unsat] = Orpp[unsat]
    Or[~unsat] = Orp[~unsat]

    unsat = D4 > 0.0

    # ------------------------------------------------------------------ #
    # Step 16: Dicalcium silicate / wollastonite
    # ------------------------------------------------------------------ #
    Cs   = _alloc()
    Wopp = _alloc()
    Wo   = _alloc()
    D5   = _alloc()

    ind = D4 < Wop / 2
    Cs[   unsat &  ind] = D4[unsat & ind]
    Wopp[ unsat &  ind] = (Wop - 2*D4)[unsat & ind]
    D5[   unsat &  ind] = 0.0

    Cs[   unsat & ~ind] = Wop[unsat & ~ind] / 2
    Wopp[ unsat & ~ind] = 0.0
    D5[   unsat & ~ind] = (D4 - Wop/2)[unsat & ~ind]

    Wo[ unsat] = Wopp[unsat]
    Wo[~unsat] = Wop[~unsat]

    unsat = D5 > 0.0

    # ------------------------------------------------------------------ #
    # Step 17: Adjust diopside / dicalcium silicate / olivine
    # ------------------------------------------------------------------ #
    Olpp = _alloc()
    Dipp = _alloc()
    Di   = _alloc()
    Ol   = _alloc()
    D6   = _alloc()

    ind = D5 < Dip
    Cs[   unsat &  ind] += (D5/2)[unsat & ind]
    Olpp[ unsat &  ind]  = (Olp + D5/2)[unsat & ind]
    Dipp[ unsat &  ind]  = (Dip - D5)[unsat & ind]
    D6[   unsat &  ind]  = 0.0

    Cs[   unsat & ~ind] += (Dip/2)[unsat & ~ind]
    Olpp[ unsat & ~ind]  = (Olp + Dip/2)[unsat & ~ind]
    Dipp[ unsat & ~ind]  = 0.0
    D6[   unsat & ~ind]  = (D5 - Dip)[unsat & ~ind]

    Ol[ unsat] = Olpp[unsat]
    Ol[~unsat] = Olp[~unsat]
    Di[ unsat] = Dipp[unsat]
    Di[~unsat] = Dip[~unsat]

    unsat = D6 > 0.0

    # ------------------------------------------------------------------ #
    # Step 18: Kaliophilite / leucite
    # ------------------------------------------------------------------ #
    Kp   = _alloc()
    Lcpp = _alloc()
    Lc   = _alloc()

    ind = Lcp >= D6 / 2
    Kp[   unsat &  ind] = (D6/2)[unsat & ind]
    Lcpp[ unsat &  ind] = (Lcp - D6/2)[unsat & ind]

    Kp[   unsat & ~ind] = Lcp[unsat & ~ind]
    Lcpp[ unsat & ~ind] = 0.0

    Lc[ unsat] = Lcpp[unsat]
    Lc[~unsat] = Lcp[~unsat]

    defsio2 = D6 - 2 * Kp   # SiO2 deficit after all adjustments

    # ------------------------------------------------------------------ #
    # Convert moles → wt% (multiply by effective molecular weight)
    # ------------------------------------------------------------------ #
    # Variable-composition minerals: weight depends on FeO/MnO split
    fe_mn_wt = _MW['FeO'] * xfeo + _MW['MnO'] * xmno   # effective Fe+Mn MW

    out = pd.DataFrame(index=data.index)
    out['cipw_quartz']          = Q   * _W['Qz']
    out['cipw_corundum']        = Cor * _W['C']
    out['cipw_orthoclase']      = Or  * _W['Or']
    out['cipw_albite']          = Ab  * _W['Ab']
    out['cipw_anorthite']       = An  * _W['An']
    out['cipw_nepheline']       = Ne  * _W['Ne']
    out['cipw_leucite']         = Lc  * _W['Lc']
    out['cipw_kaliophilite']    = Kp  * _W['Kp']
    out['cipw_acmite']          = Ac  * _W['Ac']
    out['cipw_sodium_ms']       = Ns  * _W['Ns']
    out['cipw_potassium_ms']    = Ks  * _W['Ks']
    out['cipw_diopside_mg']     = Di * xmgr * _W['Di_Mg']
    out['cipw_diopside_fe']     = Di * xfer * (fe_mn_wt + _W['Di_Fe_corr'])
    out['cipw_wollastonite']    = Wo  * _W['Wo']
    out['cipw_hypersthene_mg']  = Hy * xmgr * _W['Hy_Mg']
    out['cipw_hypersthene_fe']  = Hy * xfer * (fe_mn_wt + _W['Hy_Fe_corr'])
    out['cipw_olivine_mg']      = Ol * xmgr * _W['Ol_Mg']
    out['cipw_olivine_fe']      = Ol * xfer * (2*fe_mn_wt + _W['Ol_Fe_corr'])
    out['cipw_dicalcium_si']    = Cs  * _W['Cs']
    out['cipw_magnetite']       = Mt  * (fe_mn_wt + _W['Mt_corr'])
    out['cipw_ilmenite']        = Il  * (fe_mn_wt + _W['Il_corr'])
    out['cipw_hematite']        = Hm  * _W['He']
    out['cipw_sphene']          = Tn  * _W['Tn']
    out['cipw_perovskite']      = Pf  * _W['Pf']
    out['cipw_rutile']          = Ru  * _W['Ru']
    out['cipw_apatite']         = Ap_p * _W['Ap']
    out['cipw_cancrinite']      = Nc  * _W['Nc']
    out['cipw_calcite']         = Cc  * _W['Cc']
    out['cipw_free_p2o5']       = free_p2o5 * _MW['P2O5']
    out['cipw_free_co2']        = free_co2  * _MW['CO2']
    out['cipw_defsio2']         = -defsio2  * _MW['SiO2']   # negative = deficit

    # ------------------------------------------------------------------ #
    # Derived indices
    # ------------------------------------------------------------------ #
    salic  = out['cipw_quartz'] + out['cipw_orthoclase'] + out['cipw_albite'] + out['cipw_anorthite']
    femic  = (out['cipw_diopside_mg'] + out['cipw_diopside_fe'] +
              out['cipw_hypersthene_mg'] + out['cipw_hypersthene_fe'] +
              out['cipw_olivine_mg'] + out['cipw_olivine_fe'] +
              out['cipw_magnetite'] + out['cipw_ilmenite'] + out['cipw_hematite'])
    DI     = (out['cipw_quartz'] + out['cipw_orthoclase'] + out['cipw_albite'] +
              out['cipw_nepheline'] + out['cipw_leucite'] + out['cipw_kaliophilite'])
    An_frac = np.where(Ab + An > 0, An / (Ab + An), np.nan)

    out['cipw_salic']     = salic
    out['cipw_femic']     = femic
    out['cipw_di']        = DI               # differentiation index
    out['cipw_an_index']  = An_frac * 100    # anorthite index

    return out
