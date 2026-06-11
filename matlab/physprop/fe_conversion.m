% Estimates Fe2+ and Fe3+ from FeO_tot
function data = fe_conversion(data,fe_ratio)
% FE_CONVERSION - splits FeO_tot into FeO and Fe2O3.
%
%   data = fe_conversion(data,fe_ratio) will compute and add fields feo and
%   fe2o3 to a table data given an Fe2+/FeO_tot ratio.

data.feo = data.feo_tot .* fe_ratio;
data.fe2o3 = data.feo_tot .* (1 - fe_ratio) * 0.5 * molecularwt('Fe2O3')/molecularwt('FeO');

return