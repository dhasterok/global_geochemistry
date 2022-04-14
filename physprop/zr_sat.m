function [TZr,TZr_sd,M] = zr_sat(data);
% ZR_SAT - Computes zircon saturation temperature
%
%   [TZr,TZr_sigma,M] = zr_sat(data) compute the zircon saturation
%   temperature (in kelvins) for a granitic melt.
%
%   The following terms are required:
%       (in wt.%) 'sio2','tio2','al2o3','feo_tot','mgO','cao','na2o','k2o'
%       (in ppm) zr_ppm
%
%   A caveat about zircon thermometry. There is a danger of misuse with
%   these estimates.  One knows that zircon is saturated when inherited
%   zircons are present, but this also means that any whole rock analysis
%   will overestimate the amount of Zr in the melt.  As a result, zircon
%   saturation temperatures represent an upper bound on the temperature.
%   However, when inherited zircon is not present then the total Zr is
%   controlled by the concentration in the source, which it may have
%   exhausted and thus represents a lower bound on temperature.  Therefore
%   one must know whether inherited zircons exist or not before the utility
%   of this method may be realized.
%   pers. comm. (Martin Hand, 1 Dec 2017)
%
% Reference: Boehnke et al. (Chemical Geology, 2013)
%
% created: 30 November 2017 by D. Hasterok

oxlist = {'SiO2','TiO2','Al2O3','FeO_tot','MgO','CaO','Na2O','K2O'};

tmp = cat_mol_norm(data(:,lower(oxlist)),oxlist);
%tmp(:,oxlist) = data(:,lower(oxlist));

M = (tmp.Na + tmp.K + 2*tmp.Ca) ./ (tmp.Al .* tmp.Si);
s = 0.025;

sigma_zr = 0.1*data.zr_ppm;
sigma_M = ( ( (s*tmp.Na).^2 + (s*tmp.K).^2 + (2*s*tmp.Ca).^2 ) + ...
    (tmp.Na + tmp.K + 2*tmp.Ca).^2 .* ...
    ((s./tmp.Al).^2 + (s./tmp.Si).^2) ) ./ (tmp.Al.*tmp.Si);
sigma = [32 0.15 0.09];

ind = data.zr_ppm > 0;
data.zr_ppm(~ind) = NaN;
g = log(5e5./data.zr_ppm) + 1.16*(M - 1) + 1.48;

TZr = 10108./g; % Zr saturation temperature estimate
TZr_sd = sqrt( sigma(1)^2 + 10108^2 *( ...
    (sigma_zr./data.zr_ppm).^2 + ...
    (1.16*sigma_M).^2 + (sigma(2)*(M - 1)).^2 + ...
    sigma(3).^2 ) )./abs(g); % uncertainty

%tmp.TZr = TZr;
%tmp.TZr_sd = TZr_sd;

%writetable(tmp,'Zircon_sat_temp.csv');

return
