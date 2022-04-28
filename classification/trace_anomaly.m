function data = trace_anomaly(data,reference,varargin)
% TRACE_ANOMALY - Computes trace element anomalies
% 
%   data = trace_anomaly(data,reference) computes several common trace
%   element anomalies and adds them to the database.  The reference refers
%   to the model field in earthref.xlsx.  If there is more than one model
%   associated with reference, then add an optional layer input
%   data = trace_anomaly(data,reference,layer).
%
%   The anomalies are computed as
%
%       Nb/Nb* = Nb_N/(Th_N * La_N) where the _N denotes the element
%       concentration is normalized to the reference model.
%
%       Pb/Pb* = Pb_N/(Ce_N * Nd_N)
%
%       Eu/Eu* = Eu_N/(Sm_N * Gd_N)
%
%       Ce/Ce* = Nb_N/(Pr_N^2 / Nd_N)
%
%       Sr/Sr* = Sr_N/(Ce_N * Nd_N)

% load optional layer string
if nargin == 3
    layer = varargin{1};
else
    layer = '';
end

% load reference chemistry models
eref = readtable('earthref.xlsx');
if isempty(layer)
    iref = find(strcmpi(eref.model,reference) & eref.sigma == 0);
else
    iref = find(strcmpi(eref.model,reference) & strcmpi(eref.layer,layer) & eref.sigma == 0);
end

eref = fefix(eref);
eref.rock_type = cell([height(eref) 1]);
eref.rock_type(:) = {''};

% Nb anomaly
data.nb_anomaly = nan(height(data),1);
ind = data.nb_ppm > 0 & data.th_ppm > 0 & data.la_ppm > 0;
data.nb_anomaly(ind) = (data.nb_ppm(ind)/eref.nb_ppm(iref)) ./ ...
    sqrt( (data.th_ppm(ind)/eref.th_ppm(iref)) .* (data.la_ppm(ind)/eref.la_ppm(iref)) );

% Pb anomaly
data.pb_anomaly = nan(height(data),1);
ind = data.pb_ppm > 0 & data.ce_ppm > 0 & data.nd_ppm > 0;
data.pb_anomaly(ind) = (data.pb_ppm(ind)/eref.pb_ppm(iref)) ./ ...
    sqrt( (data.ce_ppm(ind)/eref.ce_ppm(iref)) .* (data.nd_ppm(ind)/eref.nd_ppm(iref)) );

% Eu anomaly
data.eu_anomaly = nan(height(data),1);
ind = data.eu_ppm > 0 & data.sm_ppm > 0 & data.gd_ppm > 0;
data.eu_anomaly(ind) = (data.eu_ppm(ind)/eref.eu_ppm(iref)) ./ ...
    sqrt( (data.sm_ppm(ind)/eref.sm_ppm(iref)) .* (data.gd_ppm(ind)/eref.gd_ppm(iref)) );

%data.eu_anomaly2 = nan(height(data),1);
%ind = data.eu_ppm > 0 & data.sm_ppm > 0 & data.gd_ppm > 0;
%data.eu_anomaly2(ind) = 2*(data.eu_ppm(ind)/eref.eu_ppm(iref)) ./ ...
%    ( (data.sm_ppm(ind)/eref.sm_ppm(iref)) + (data.gd_ppm(ind)/eref.gd_ppm(iref)) );

% Ce anomaly
% Lawrence et al. (Aquat. Geochem. 2006) 10.1007/s10498-005-4471-8
% this definition as opposed to the one below is less susceptible to La
% enrichment anomalies that can skew the Ce anomaly
data.ce_anomaly = nan(height(data),1);
ind = data.ce_ppm > 0 & data.pr_ppm > 0 & data.nd_ppm > 0;
data.ce_anomaly(ind) = (data.ce_ppm(ind)/eref.ce_ppm(iref)) ./ ...
    ( (data.pr_ppm(ind)/eref.pr_ppm(iref)).^2 ./ (data.nd_ppm(ind)/eref.nd_ppm(iref)) );

%data.ce_anomaly2 = nan(height(data),1);
%ind = data.ce_ppm > 0 & data.pr_ppm > 0 & data.la_ppm > 0;
%data.ce_anomaly2(ind) = 2*(data.ce_ppm(ind)/eref.ce_ppm(iref)) ./ ...
%    ( (data.la_ppm(ind)/eref.la_ppm(iref)) + (data.pr_ppm(ind)/eref.pr_ppm(iref)) );

% Sr anomaly
data.sr_anomaly = nan(height(data),1);
ind = data.sr_ppm > 0 & data.ce_ppm > 0 & data.nd_ppm > 0;
data.sr_anomaly(ind) = (data.sr_ppm(ind)/eref.sr_ppm(iref)) ./ ...
    sqrt( (data.ce_ppm(ind)/eref.ce_ppm(iref)) .* (data.nd_ppm(ind)/eref.nd_ppm(iref)) );


return