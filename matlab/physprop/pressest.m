function [P_qtz,P_fsp] = pressest(data)
% Pressure estimate of haplogranite system from Yang [Lithos, 2017]

% step 1 - normalize to anhydrous conditions
oxides = {'SiO2'; 'TiO2'; 'Al2O3'; 'FeO'; 'MnO'; ...
    'MgO'; 'CaO'; 'Na2O'; 'K2O'; 'P2O5'};
tmp = oxide_norm(data,oxides);

cipw = cipwnorm(data);

% calibrated range
ind = 15 <= cipw.Quartz & cipw.Quartz <= 40;

T = cipw.Quartz + cipw.Albite + cipw.Orthoclase;
Q = 100*cipw.Quartz./T;
Ab = 100*cipw.Albite./T;
Or = 100*cipw.Orthoclase./T;

P_qtz = -0.2426*Q.^3 + 26.392*Q.^2 - 980.74*Q + 12563;
P_qtz(~ind) = NaN;


x = Ab + Or;
% calibrated range
ind = 60 <= x & x <= 85;

P_fsp = 0.2426*x.^3 - 46.397*x.^2 + 2981.3*x - 64224;
P_fsp(~ind) = NaN;

return
