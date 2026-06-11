function M = mgnum(MgO,FeO,Fe2O3)
% MGNUM - Magnesium number
%
% Computes the magnesium number,
%                 MgO
% Mg# = 100 x -----------
%             (MgO + FeO)
% where MgO and FeO are in moles.  MGNUM accepts the mass fraction of each.  For
% rocks with FeOT, 0.1% is assumed to be Fe2O3.

Fe = nan(size(MgO));
ind = ((Fe2O3 == 0 | isnan(Fe2O3)) & FeO ~= 0);
Fe(ind) = FeO(ind)/molecularwt('FeO');

ind = (Fe2O3 ~= 0 & (isnan(FeO) | FeO == 0));
Fe(ind) = 2*0.89*Fe2O3(ind)/molecularwt('Fe2O3');

ind = (Fe2O3 ~= 0 & FeO ~= 0);
if sum(ind) > 0
    Fe(ind) =  0.89*FeO/molecularwt('FeO');
end

Mg = MgO/molecularwt('MgO');

M = Mg./(Mg + Fe)*100;


return
