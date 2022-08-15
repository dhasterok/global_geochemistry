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
data.thermal_conductivity = exp(4.899*data.sio2 + 0.01736*data.tio2 + ...
    0.8895*data.al2o3 + 6.641*data.feo_tot + 2.136*data.mgo + 0.7416*data.cao + ...
    0.05848*data.na2o + 0.2804*data.k2o);

ind = (data.thermal_conductivity < 1.5 | 5 < data.thermal_conductivity);
data.thermal_conductivity(ind) = nan;

return