function data = tcest(data)
% TCEST - estimates thermal conductivity from a data table
%
% data = tcest(data) computes the heat production for a table adding
% a column for thermal conductivity (W m^-1 K^-1)
%
% results from Jennings et al., GJI, (2019)
%oxides = {'SiO2'; 'TiO2'; 'Al2O3'; 'FeO'; 'MgO'; 'CaO'; 'Na2O'; 'K2O'; 'LOI'};
%tmp = oxide_norm(data,oxides);

%data.thermal_conductivity = exp( (1.729*tmp.sio2 - 0.253*tmp.al2o3 + ...
%    1.054*tmp.feo_tot - 3.831*tmp.na2o - 1.622*tmp.k2o + ...
%    1.737*tmp.loi) /100 );

%data.thermal_conductivity = exp( (1.72*tmp.sio2 + 1.018*tmp.mgo ...
%    - 3.652*tmp.na2o - 1.791*tmp.k2o)/100 );

% results from Willcocks et al. (in prep)
oxides = {'SiO2'; 'TiO2'; 'Al2O3'; 'FeO'; 'MgO'; 'CaO'; 'Na2O'; 'K2O'};
oxlist = {'sio2'; 'tio2'; 'al2o3'; 'feo_tot'; 'mgo'; 'cao'; 'na2o'; 'k2o'};
tmp = oxide_norm(data,'Normalization','anhydrous','Oxides',oxlist);

data.thermal_conductivity = exp((1.5542*tmp.sio2 - 6.5211*tmp.tio2 ...
    - 0.2547*tmp.al2o3 + 2.0660*tmp.feo_tot + 0.5472*tmp.mgo ...
    + 0.2258*tmp.cao - 2.4017*tmp.na2o - 0.6721*tmp.k2o)./sum(tmp{:,oxlist},2));

ind = (data.thermal_conductivity < 1.5 | 5 < data.thermal_conductivity);
data.thermal_conductivity(ind) = nan;

return