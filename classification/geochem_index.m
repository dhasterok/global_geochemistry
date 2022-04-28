function data = geochem_index(data)
% GEOCHEM_INDEX - Computes several popular geochemical classification
% indicies (specifically for igneous rocks)
%
% Assumes all iron is reported as FeO(T).
% Frost et al. 2001

% initialize indicies as NaN
data.Mg_number = nan([height(data) 1]);
data.Fe_number = nan([height(data) 1]);
data.MALI = nan([height(data) 1]);
data.ASI = nan([height(data) 1]);
data.maficity = nan([height(data) 1]);
data.CIA = nan([height(data) 1]);
data.WIP = nan([height(data) 1]);
data.CPA = nan([height(data) 1]);
data.spar = nan([height(data) 1]);
data.qtzindex = nan([height(data) 1]);
data.R1 = nan([height(data) 1]);
data.R2 = nan([height(data) 1]);

% convert all BDL to 0 (assume they are small enough that they won't affect
% ratios
sio2 = data.sio2;
sio2(sio2 < 0) = 0;
tio2 = data.tio2;
tio2(tio2 < 0) = 0;
tio2(isnan(tio2)) = 0;
al2o3 = data.al2o3;
al2o3(al2o3 < 0) = 0;
feo = data.feo_tot;
feo(feo < 0) = 0;
mgo = data.mgo;
mgo(mgo < 0) = 0;
cao = data.cao;
cao(cao < 0) = 0;
na2o = data.na2o;
na2o(na2o < 0) = 0;
k2o = data.k2o;
k2o(k2o < 0) = 0;
p2o5 = data.p2o5;
p2o5(p2o5 < 0) = 0;
p2o5(isnan(p2o5)) = 0;

% mole fraction
n_sio2 = sio2/molecularwt('SiO2');
n_tio2 = tio2/molecularwt('TiO2');
n_al2o3 = al2o3/molecularwt('Al2O3');
n_feo = feo/molecularwt('FeO');
n_mgo = mgo/molecularwt('MgO');
n_cao = cao/molecularwt('CaO');
n_na2o = na2o/molecularwt('Na2O');
n_k2o = k2o/molecularwt('K2O');
n_p2o5 = p2o5/molecularwt('P2O5');

% reduce CaO for apatite
n_canoap = n_cao - 10/3*n_p2o5;
n_canoap(n_canoap < 0) = 0;

% Mg number
% -------------------------------------
% computation of fe2+ from other oxides
% Le Maitre, R.W. (1976) The Chemical Variability of Some Common Igneous
% Rocks. Journal of Petrology, 17, 589--598,
% doi:10.1093/petrology/17.4.589.
data.fe2_fe_est = nan([height(data) 1]);
ind = rockgroup(data,'all plutonic');
if sum(ind) > 0
    data.fe2_fe_est(ind) = 0.93 - 0.0042*sio2(ind) - 0.022*(na2o(ind) + k2o(ind));
end
ind = rockgroup(data,'all volcanic');
if sum(ind) > 0
    data.fe2_fe_est(ind) = 0.88 - 0.0016*sio2(ind) - 0.027*(na2o(ind) + k2o(ind));
end


data.Mg_number = n_mgo ./ (n_mgo + 0.85*n_feo);


% Fe number
% -------------------------------------
% Frost et al. (2001)
data.Fe_number = feo ./ (feo + mgo);
data.Fe_number(isinf(data.Fe_number)) = NaN;

% (MALI) Modified Alkali Lime Index
% -------------------------------------
% Frost et al. (2001)

data.MALI = na2o + k2o - cao;

% (CAI) Calcic-Alkalic Index
% -------------------------------------
% Sets as zero the data along the calc-alkalic/alkali-calcic dividing line.
% Calcic is approximately -? and alkalic +?.
data.CAI = data.MALI - (-44.72 + 1.094*sio2 - 0.00527*sio2.^2);

% ASI - Alumina Saturation Index
% -------------------------------------
% Frost et al. (2001)
data.ASI = n_al2o3 ./ ( n_canoap + (n_na2o + n_k2o) );
data.ASI(isinf(data.ASI)) = NaN;

% AI - Agpaitic Index
% -------------------------------------
% 
data.AI = n_al2o3 ./ ( n_na2o + n_k2o );
data.AI(isinf(data.AI)) = NaN;

% Maficity
% -------------------------------------
% From Debon & LeFort (Trans. R. Soc. Edinb. Earth Sci., 1983) and
% Debon and LeFort (Bull. Mineral, 1988)
data.maficity = n_mgo + n_feo + n_tio2;


% Feldspar index
% -------------------------------------
% From Debon & LeFort (Trans. R. Soc. Edinb. Earth Sci., 1983) and
% Debon and LeFort (Bull. Mineral, 1988)
% Modified to remove apatite
data.spar = 2*n_k2o - (2*n_na2o + n_cao - 10/3*n_p2o5);


% Quartz index
% -------------------------------------
% From Debon & LeFort (Trans. R. Soc. Edinb. Earth Sci., 1983) and
% Debon and LeFort (Bull. Mineral, 1988)
data.qtzindex = n_sio2/3 - (2*n_na2o + 2*n_k2o + (n_canoap)/3);


% CIA - Chemical index of alteration
% -------------------------------------
% from Price and Velbel (Chem. Geol., 2003)
data.CIA = 100 * n_al2o3 ./ ( n_al2o3 + n_canoap + n_na2o + n_k2o );
data.CIA(isinf(data.CIA)) = NaN;

% Parker weathering data (allows for Al mobility)
% -------------------------------------
% from Price and Velbel (Chem. Geol., 2003)
% optimum fresh value > 100
% optimum weathered value 0
% change to molar proportions?
data.WIP = 100*(2*n_na2o/0.35 + n_mgo/0.9 + 8*n_k2o + n_cao/0.7);


% CPA - Chemical proxy of alteration
% -------------------------------------
% Buggle (Quaternary Int., 2011) doi: 10.1016/j.quaint.2010.07.019
data.CPA = 100*n_al2o3 ./ ( n_al2o3 + n_na2o );
data.CPA(isinf(data.CPA)) = NaN;

% R1 & R2
% ------------------------------------
% De La Roche et al. (1980)
data.R1 = 4*n_sio2 - 22*(n_na2o + n_k2o) - 2*(n_feo + n_tio2);
data.R2 = 6*n_cao + 2*n_mgo + 2*n_al2o3;

return
