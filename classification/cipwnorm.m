function [cipw,P] = cipwnorm(data, varargin)
% CIPW norm (Verma, 2002, 2003)
%------------------------------------

% default parameters
volc = logical(ones([height(data) 1]));
AdjTAS = 0;
Cancrinite = 0;
Calcite = 0;

if nargin > 2
    volc = logical(varargin{1});
elseif any(strcmp(data.Properties.VariableNames,'rock_group'))
    volc = rockgroup(data,'all volcanic');
end

if ~any(strcmp(data.Properties.VariableNames,'co2'));
    data.co2 = zeros([size(data.sio2)]);
end

% Calculate major elements data on an anhydrous basis
%------------------------------------
data = adj_rock(data, volc, AdjTAS, Cancrinite, Calcite);


% Load chemical data
%------------------------------------
MinWeight = readtable('minweight.csv','Format','%s%s%s%f%f%f');
MinWeight.Properties.RowNames = MinWeight.MineralName;
fmt = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f';
OxiWeight = readtable('oxiweight.csv','Format',fmt);

oxlist = {'sio2','tio2','al2o3','mno','feo','fe2o3','mgo','cao', ...
    'na2o','k2o','p2o5','co2','Fe_ratio'};

tmp = data(:,oxlist);

% Mole computation
%------------------------------------
tmp.n_sio2 = tmp.sio2/molecularwt('SiO2');
tmp.n_tio2 = tmp.tio2/molecularwt('TiO2');
tmp.n_al2o3 = tmp.al2o3/molecularwt('Al2O3');
tmp.n_mno = data.mno/molecularwt('MnO'); 
tmp.n_feo = tmp.feo/molecularwt('FeO') + tmp.mno/molecularwt('MnO');
tmp.n_fe2o3 = tmp.fe2o3/molecularwt('Fe2O3');
tmp.n_mgo = tmp.mgo/molecularwt('MgO');
tmp.n_cao = tmp.cao/molecularwt('CaO');
tmp.n_na2o = tmp.na2o/molecularwt('Na2O'); 
tmp.n_k2o = tmp.k2o/molecularwt('K2O');
tmp.n_p2o5 = tmp.p2o5/molecularwt('P2O5'); 
tmp.xmno = tmp.n_mno./tmp.n_feo;
tmp.xfeo = (tmp.feo/molecularwt('FeO'))./tmp.n_feo;
tmp.n_co2 = tmp.co2/molecularwt('CO2'); 

cipw = tmp;
% Main functions
%------------------------------------
normlist = {'Ap','FREE_P2O5','Nc','Cc','FREE_CO2','Il','Orp','Ks','Y', ...
    'Abp','Ac','Ns','An','C','Tnp','Ru','Mt','Hm','Dip','Wop', ...
    'Hyp','Q','D','Olp','Hy','D1','Pf','Tn','D2','Ne','Abpp','D3', ...
    'Ab','Lcp','Orpp','D4','Or','Cs','Wopp','D5','Wo','Olpp', ...
    'Dipp','D6','Ol','Di','Kp','Lcpp','Lc','DEFSIO2'};

m = array2table(zeros([height(data) length(normlist)]), ...
    'VariableNames',normlist);

% Apatite/FREE_P2O5
%------------------------------------
ind = tmp.n_cao >= 10/3 * tmp.n_p2o5;
m.Ap(ind) = tmp.n_p2o5(ind);
tmp.n_cao(ind) = tmp.n_cao(ind) - 10/3 * m.Ap(ind);
tmp.n_p2o5(ind) = 0;

m.Ap(~ind) = 3/10 * tmp.n_cao(~ind);
tmp.n_p2o5(~ind) = tmp.n_p2o5(~ind) - m.Ap(~ind);
tmp.n_cao(~ind) = 0;

m.p2o5 = tmp.n_p2o5;

% Sodium carbonate/calcite
%------------------------------------
m.Nc(~Cancrinite) = 0;

ind = tmp.n_na2o >= tmp.n_co2;
m.Nc(Cancrinite & ind) = tmp.n_co2(Cancrinite & ind);
tmp.n_na2o(Cancrinite & ind) = tmp.n_na2o(Cancrinite & ind) - ...
    tmp.n_co2(Cancrinite & ind);
tmp.n_co2(Cancrinite & ind) = 0;

m.Nc(Cancrinite & ~ind) = tmp.n_na2o(Cancrinite & ~ind);
tmp.n_na2o(Cancrinite & ~ind) = 0;
tmp.n_co2(Cancrinite & ~ind) = tmp.n_co2(Cancrinite & ~ind) - ...
    m.Nc(Cancrinite & ~ind);
tmp.n_na2o(Cancrinite & ~ind) = 0;

% Calcite
%------------------------------------
m.Cc(~Calcite) = 0;

ind = tmp.n_cao >= tmp.n_co2;
m.Cc(Calcite & ind) = tmp.n_co2(Calcite & ind);
tmp.n_cao(Calcite & ind) = tmp.n_cao(Calcite & ind) - ...
    m.Cc(Calcite & ind);
tmp.n_co2(Calcite & ind) = 0;

m.Cc(Calcite & ~ind) = tmp.n_cao(Calcite & ~ind);
tmp.n_cao(Calcite & ~ind) = 0;
tmp.n_co2(Calcite & ~ind) = tmp.n_co2(Calcite & ~ind) - ...
    m.Cc(Calcite & ~ind);
tmp.n_cao(Calcite & ~ind) = 0;

m.co2 = tmp.n_co2;

% Ilmenite
%------------------------------------
ind = tmp.n_feo >= tmp.n_tio2;
m.Il(ind) = tmp.n_tio2(ind);
tmp.n_feo(ind) = tmp.n_feo(ind) - m.Il(ind);
tmp.n_tio2(ind) = 0;

m.Il(~ind) = tmp.n_feo(~ind);
tmp.n_tio2(~ind) = tmp.n_tio2(~ind) - m.Il(~ind);
tmp.n_feo(~ind) = 0;

% Orthoclase/Potassium metasilicate
%------------------------------------
ind = tmp.n_al2o3 >= tmp.n_k2o;
m.Orp(ind) = tmp.n_k2o(ind);
tmp.n_al2o3(ind) = tmp.n_al2o3(ind) - m.Orp(ind);
tmp.n_k2o(ind) = 0;

m.Orp(~ind) = tmp.n_al2o3(~ind);
tmp.n_k2o(~ind) = tmp.n_k2o(~ind) - m.Orp(~ind);
tmp.n_al2o3(~ind) = 0;

m.Ks = tmp.n_k2o;
m.Y = 6 * m.Orp + m.Ks;

% Albite
%------------------------------------ 
ind = tmp.n_al2o3 >= tmp.n_na2o;
m.Abp(ind) = tmp.n_na2o(ind);
tmp.n_al2o3(ind) = tmp.n_al2o3(ind) - m.Abp(ind);
tmp.n_na2o(ind) = 0;

m.Abp(~ind) = tmp.n_al2o3(~ind);
tmp.n_na2o(~ind) = tmp.n_na2o(~ind) - m.Abp(~ind);
tmp.n_al2o3(~ind) = 0;

m.Y = m.Y + 6 * m.Abp;

% Acmite/Sodium metasilicate
%------------------------------------ 
ind = tmp.n_na2o >= tmp.n_fe2o3;
m.Ac(ind) = tmp.n_fe2o3(ind);
tmp.n_na2o(ind) = tmp.n_na2o(ind) - m.Ac(ind);
tmp.n_fe2o3(ind) = 0;

m.Ac(~ind) = tmp.n_na2o(~ind);
tmp.n_fe2o3(~ind) = tmp.n_fe2o3(~ind) - m.Ac(~ind);
tmp.n_na2o(~ind) = 0;

m.Ns = tmp.n_na2o;
m.Y = 4 * m.Ac + m.Ns + m.Y;

% Anorthite/Corundum
%------------------------------------ 
ind = tmp.n_al2o3 >= tmp.n_cao;
m.An(ind) = tmp.n_cao(ind);
tmp.n_al2o3(ind) = tmp.n_al2o3(ind) - m.An(ind);
tmp.n_cao(ind) = 0;
    
m.An(~ind) = tmp.n_al2o3(~ind);
tmp.n_cao(~ind) = tmp.n_cao(~ind) - m.An(~ind);
tmp.n_al2o3(~ind) = 0;

m.C = tmp.n_al2o3;
m.Y = 2 * m.An + m.Y;


% Sphene/Rutile
%------------------------------------ 
ind = tmp.n_cao >= tmp.n_tio2;
m.Tnp(ind) = tmp.n_tio2(ind);
tmp.n_cao(ind) = tmp.n_cao(ind) - m.Tn(ind);
tmp.n_tio2(ind) = 0;

m.Tnp(~ind) = tmp.n_cao(~ind);
tmp.n_tio2(~ind) = tmp.n_tio2(~ind) - m.Tnp(~ind);
tmp.n_cao(~ind) = 0;

m.Ru = tmp.n_tio2;
m.Y = m.Tnp + m.Y;

% Magnetite/Hematite
%------------------------------------
ind = tmp.n_fe2o3 >= tmp.n_feo;
m.Mt(ind) = tmp.n_feo(ind);
tmp.n_fe2o3(ind) = tmp.n_fe2o3(ind) - m.Mt(ind);
tmp.n_feo(ind) = 0;

m.Mt(~ind) = tmp.n_fe2o3(~ind);
tmp.n_feo(~ind) = tmp.n_feo(~ind) - m.Mt(~ind);
tmp.n_fe2o3(~ind) = 0;

m.Hm = tmp.n_fe2o3;

% Subdivision of some normative minerals (end-members)
%------------------------------------%
tmp.n_femg = tmp.n_feo + tmp.n_mgo;

tmp.xfer = tmp.n_feo ./ (tmp.n_feo + tmp.n_mgo);
tmp.xmgr = tmp.n_mgo ./ (tmp.n_mgo + tmp.n_feo);

% Diopside/Wollastonite/Hypersthene
%------------------------------------ 
ind = tmp.n_cao >= tmp.n_femg;
m.Dip(ind) = tmp.n_femg(ind);
tmp.n_cao(ind) = tmp.n_cao(ind) - m.Dip(ind);
tmp.n_femg(ind) = 0;

m.Dip(~ind) = tmp.n_cao(~ind);
tmp.n_femg(~ind) = tmp.n_femg(~ind) - m.Dip(~ind);
tmp.n_cao(~ind) = 0;

m.Wop = tmp.n_cao;
m.Hyp = tmp.n_femg;

ind = m.Wop > 0;
m.Y = 2 * m.Dip + m.Wop + m.Hyp + m.Y;

% Quartz/Undersatured Quartz
%------------------------------------ 
ind = tmp.n_sio2 >= m.Y;
m.Q(ind) = tmp.n_sio2(ind) - m.Y(ind);
m.D(ind) = 0;

m.D(~ind) = m.Y(~ind) - tmp.n_sio2(~ind);
m.Q(~ind) = 0;

unsaturated = ~ind;

% Olivine/Hypersthene
%------------------------------------ 
ind = m.D < m.Hyp/2;
m.Olp(ind) = m.D(ind);
m.Hy(ind) = m.Hyp(ind) - 2 * m.D(ind);
m.D1(ind) = 0;

m.Olp(~ind) = m.Hyp(~ind)/2;
m.Hy(~ind) = 0;
m.D1(~ind) = m.D(~ind) - m.Hyp(~ind)/2;

m.Olp(~unsaturated) = 0;
m.Hy(~unsaturated) = m.Hyp(~unsaturated);
m.D1(~unsaturated) = 0;

unsaturated = m.D1 > 0;

% Sphene/Perovskite
%------------------------------------ 
ind = m.D1 < m.Tnp;
m.Pf(ind) = m.D1(ind);
m.Tn(ind) = m.Tnp(ind) - m.D1(ind);
m.D2(ind) = 0;

m.Pf(~ind) = m.Tnp(~ind);
m.Tn(~ind) = 0;
m.D2(~ind) = m.D1(~ind) - m.Tnp(~ind);

m.Pf(~unsaturated) = 0;

unsaturated = m.D2 > 0;

% Nepheline/Albite
%------------------------------------ 
ind = m.D2 < 4*m.Abp;
m.Ne(ind) = m.D2(ind)/4;
m.Abpp(ind) = m.Abp(ind) - m.D2(ind)/4;
m.D3(ind) = 0;

m.Ne(~ind) = m.Abp(~ind);
m.Abpp(~ind) = 0;
m.D3(~ind) = m.D2(~ind) - m.Abp(~ind)*4;

m.Ne(~unsaturated) = 0;
m.Ab(unsaturated) = m.Abpp(unsaturated);
m.Ab(~unsaturated) = m.Abp(~unsaturated);

unsaturated = m.D3 > 0;


% Leucite/Orthoclase
%------------------------------------ 
ind = m.D3 < 2*m.Orp;
m.Lcp(ind) = m.D3(ind)/2;
m.Orpp(ind) = m.Orp(ind) - m.D3(ind)/2;
m.D4(ind) = 0;

m.Lcp(~ind) = m.Orp(~ind);
m.Orpp(~ind) = 0;
m.D4(~ind) = m.D3(~ind) - 2 * m.Orp(~ind);

m.Lcp(~unsaturated) = 0;
m.Or(unsaturated) = m.Orpp(unsaturated);
m.Or(~unsaturated) = m.Orp(~unsaturated);

unsaturated = m.D4 > 0;

% Dicalcium silicate/Wollastonite
%------------------------------------
ind = m.D4 < m.Wop/2;
m.Cs(ind) = m.D4(ind);
m.Wopp(ind) = m.Wop(ind) - 2 * m.D4(ind);
m.D5(ind) = 0;

m.Cs(~ind) = m.Wop(~ind)/2;
m.Wopp(~ind) = 0;
m.D5(~ind) = m.D4(~ind) - m.Wop(~ind)/2;

m.Cs(~unsaturated) = 0;
m.Wo(unsaturated) = m.Wopp(unsaturated);
m.Wo(~unsaturated) = m.Wop(~unsaturated);

unsaturated =  m.D5 > 0;

% Adjust Diopside/Dicalcium silicate/Olivine
%------------------------------------ 
ind = m.D5 < m.Dip;
m.Cs(ind) = m.Cs(ind) + m.D5(ind)/2;
m.Olpp(ind) = m.Olp(ind) + m.D5(ind)/2;
m.Dipp(ind) = m.Dip(ind) - m.D5(ind);
m.D6(ind) = 0;

m.Cs(~ind) = m.Cs(~ind) + m.Dip(~ind)/2;
m.Olpp(~ind) = m.Olp(~ind) + m.Dip(~ind)/2;
m.Dipp(~ind) = 0;
m.D6(~ind) = m.D5(~ind) - m.Dip(~ind);

m.Ol(unsaturated) = m.Olpp(unsaturated);
m.Ol(~unsaturated) = m.Olp(~unsaturated);
m.Di(unsaturated) = m.Dipp(unsaturated);
m.Di(~unsaturated) = m.Dip(~unsaturated);

unsaturated = m.D6 > 0;


% Adjust Kaliophilite/Leucite
%------------------------------------ 
ind = m.Lcp >= m.D6/2;
m.Kp(ind) = m.D6(ind)/2;
m.Lcpp(ind) = m.Lcp(ind) - m.D6(ind)/2;

m.Kp(~ind) = m.Lcp(~ind);
m.Lcpp(~ind) = 0;

m.Lc(unsaturated) = m.Lcpp(unsaturated);
m.Lc(~unsaturated) = m.Lcp(~unsaturated);
m.DEFSIO2  = m.D6 - 2 * m.Kp;

% Print Minerals
%------------------------------------%


cipw.Quartz = m.Q * MinWeight{'Quartz','ConsWeight'};
cipw.Corundum = m.C * MinWeight{'Corundum','ConsWeight'};
cipw.Orthoclase = m.Or * MinWeight{'Orthoclase','ConsWeight'};
cipw.Albite = m.Ab * MinWeight{'Albite','ConsWeight'};
cipw.Anorthite = m.An * MinWeight{'Anorthite','ConsWeight'};
cipw.Nepheline = m.Ne * MinWeight{'Nepheline','ConsWeight'};
cipw.Leucite = m.Lc * MinWeight{'Leucite','ConsWeight'};
cipw.Kaliophilite = m.Kp * MinWeight{'Kaliophilite','ConsWeight'};
cipw.Acmite = m.Ac * MinWeight{'Acmite','ConsWeight'};
cipw.Sodium_metasilicate = m.Ns * MinWeight{'Sodium_metasilicate','ConsWeight'};
cipw.Potassium_metasilicate = m.Ks * MinWeight{'Potassium_metasilicate','ConsWeight'};
cipw.Diopside_Mg = m.Di .* tmp.xmgr * MinWeight{'Diopside_Mg','ConsWeight'};
cipw.Diopside_Fe = m.Di .* tmp.xfer .* ...
    (molecularwt('FeO') * tmp.xfeo + molecularwt('MnO')*tmp.xmno + ...
    MinWeight{'Diopside_Fe','CorrWeight'});
cipw.Wollastonite = m.Wo * MinWeight{'Wollastonite','ConsWeight'};
cipw.Hypersthene_Mg = m.Hy .* tmp.xmgr * MinWeight{'Hypersthene_Mg','ConsWeight'};
cipw.Hypersthene_Fe = m.Hy .* tmp.xfer .* ...
    (molecularwt('FeO') * tmp.xfeo + ...
    molecularwt('MnO') * tmp.xmno +...
    MinWeight{'Hypersthene_Fe','CorrWeight'});
cipw.Olivine_Mg = m.Ol .* tmp.xmgr * MinWeight{'Olivine_Mg','ConsWeight'};
cipw.Olivine_Fe = m.Ol .* tmp.xfer .* ...
    (molecularwt('FeO') * tmp.xfeo*2 + ...
    molecularwt('MnO') * tmp.xmno*2 + ...
    MinWeight{'Olivine_Fe','CorrWeight'});
cipw.Dicalcium_silicate = m.Cs * MinWeight{'Dicalcium_silicate','ConsWeight'};
cipw.Magnetite = m.Mt .* ...
    (molecularwt('FeO') * tmp.xfeo + ...
    molecularwt('MnO') * tmp.xmno + ...
    MinWeight{'Magnetite','CorrWeight'});
cipw.Ilmenite = m.Il .* ...
    (molecularwt('FeO') * tmp.xfeo + ...
    molecularwt('MnO') * tmp.xmno + ...
    MinWeight{'Ilmenite','CorrWeight'});
cipw.Hematite = m.Hm * MinWeight{'Hematite','ConsWeight'};
cipw.Sphene = m.Tn * MinWeight{'Sphene','ConsWeight'};
cipw.Perovskite = m.Pf * MinWeight{'Perovskite','ConsWeight'};
cipw.Rutile = m.Ru * MinWeight{'Rutile','ConsWeight'};
cipw.Apatite_Ca = m.Ap * MinWeight{'Apatite_Ca','ConsWeight'};
cipw.Sodium_Carbonate = m.Nc * MinWeight{'Sodium_Carbonate','ConsWeight'};
cipw.Calcite = m.Cc * MinWeight{'Calcite','ConsWeight'};
cipw.DEFSIO2 = -m.DEFSIO2 * molecularwt('SiO2');
cipw.FREE_P2O5 = m.p2o5 * molecularwt('P2O5');
cipw.FREE_CO2 = m.co2 * molecularwt('CO2');
cipw.An_index = m.An ./ (m.Ab + m.An)*100;

% Geochemical Indicies
%------------------------------------
P.Salic = cipw.Quartz + cipw.Orthoclase + ...
    cipw.Albite + cipw.Anorthite;
P.Femic = cipw.Diopside_Mg + cipw.Diopside_Fe + ...
    cipw.Hypersthene_Mg + cipw.Hypersthene_Fe + ...
    cipw.Olivine_Mg + cipw.Olivine_Fe + ...
    cipw.Magnetite + cipw.Ilmenite + ...
    cipw.Hematite;
P.CI = cipw.Anorthite + (2.1570577*cipw.Diopside_Mg) + ...
    cipw.Olivine_Mg + (0.7007616*cipw.Hypersthene_Fe);
P.DI =  cipw.Quartz + cipw.Orthoclase + ...
    cipw.Albite + cipw.Nepheline + cipw.Leucite;
P.SI = (100*tmp.n_mgo)./(tmp.n_mgo + tmp.n_fe2o3 + ...
    tmp.n_feo + tmp.n_na2o + tmp.n_k2o);
ind =  tmp.n_sio2 > 0 & tmp.n_k2o./tmp.n_na2o > 1.0 & ...
    tmp.n_k2o./tmp.n_na2o < 2.5;
P.AR = (tmp.n_al2o3(ind) + tmp.n_cao(ind) + (2*tmp.n_na2o(ind))) ./ ...
    (tmp.n_al2o3(ind) + tmp.n_cao(ind) - 2*tmp.n_na2o(ind));
    P.AR = (tmp.n_al2o3(~ind) + tmp.n_cao(~ind) + tmp.n_na2o(~ind) + tmp.n_k2o(~ind))./ ...
    (tmp.n_al2o3(~ind) + tmp.n_cao(~ind) - tmp.n_na2o(~ind) - tmp.n_k2o(~ind));

P.Mg_number = 100 * ((molecularwt('Mg')/molecularwt('MgO')) * tmp.n_mgo)./ ...
    (((molecularwt('Mg')/molecularwt('MgO')) * tmp.n_mgo) + ...
    ((molecularwt('Fe')/molecularwt('FeO')) * tmp.n_feo));
P.Fe_number = tmp.n_feo ./ (tmp.n_feo + tmp.n_mgo);
P.MALI =  tmp.n_na2o + tmp.n_k2o - tmp.n_cao;
P.ACNK = ((molecularwt('Al')/molecularwt('Al2O3')) * tmp.n_al2o3) ./ ...
    ((molecularwt('Ca')/molecularwt('CaO') * tmp.n_cao - ...
    1.67 * molecularwt('P') / molecularwt('P2O5') * tmp.n_p2o5) + ...
    molecularwt('Na') / molecularwt('Na2O') * tmp.n_na2o + ...
    molecularwt('K') / molecularwt('K2O') * tmp.n_k2o);
P.ANK = molecularwt('Al')/molecularwt('Al2O3') * tmp.n_al2o3 ./ ...
    (molecularwt('Na')/molecularwt('Na2O') * tmp.n_na2o + ...
    molecularwt('K')/molecularwt('K2O') * tmp.n_k2o);
P.AI = molecularwt('Al') / molecularwt('Al2O3') * tmp.n_al2o3 - ...
    molecularwt('K') / molecularwt('K2O') * tmp.n_k2o - ...
    molecularwt('Na')/molecularwt('Na2O') * tmp.n_na2o;
P.FSSI = cipw.Quartz - (cipw.Leucite  + ...
    2*cipw.Nepheline + cipw.Kaliophilite) /100;


% Calculate the theoretical density of rocks
%------------------------------------

function density = rockdensity(cipw,MinWeight)

% Matrix computation
%------------------------------------

minlist = {'Quartz',
    'Corundum',
    'Orthoclase',
    'Albite',
    'Anorthite',
    'Nepheline',
    'Leucite',
    'Kaliophilite',
    'Acmite',
    'Sodium_metasilicate',
    'Potassium_metasilicate',
    'Apatite',
    'Diopside',
    'Hedenbergite',
    'Wollastonite',
    'Enstatite',
    'Ferrosilite',
    'Forsterite',
    'Fayalite',
    'Dicalcium_silicate',
    'Magnetite',
    'Ilmenite',
    'Hematite',
    'Sphene',
    'Perovskite',
    'Rutile',
    'Sodium_Carbonate',
    'Calcite'};

% Vol. prop. of mineral
%------------------------------------
cipw.volume(minlist) = cipw.norm(minlist) ./ MinWeight(minlist,'Density');

total_volume = sum(cipw.volume(minlist));

% Density of minerals
%------------------------------------
cipw.density(minlist) = cipw.volume(minlist) * MinWeight(minlist,'Density') / ...
    total_volume;

Density.Rock = sum(Dens(:));

return
