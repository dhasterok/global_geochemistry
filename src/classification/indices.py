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
