function data = fefix(data);
% FEFIX - computes total iron and fe2+ to total iron ratio.
%
%   data = fefix(data) converts all iron to total iron and computes the
%   Fe2+/(Fe2+ + Fe3+) ratio.  All other Fe columns are deleted (FeO,
%   Fe2O3, and Fe2O3_tot).
%
%   To restore FeO and Fe2O3:
%
%       FeO = molecularwt('FeO')*data.fe2_fe_tot.*data.feo_tot
%       Fe2O3 = 0.5*molecularwt('Fe2O3')*(1 - data.fe2_fe_tot).*data.feo_tot
%

% conversion factor for Fe2O3 to FeO
fe_factor = 2*molecularwt('FeO')/molecularwt('Fe2O3');

if ~any(strcmp(data.Properties.VariableNames,'feo')) & ...
    ~any(strcmp(data.Properties.VariableNames,'fe2o3')) & ...
    ~any(strcmp(data.Properties.VariableNames,'fe2o3_tot'))
    return
elseif any(strcmp(data.Properties.VariableNames,'feo')) & ...
    ~any(strcmp(data.Properties.VariableNames,'fe2o3')) & ...
    ~any(strcmp(data.Properties.VariableNames,'fe2o3_tot'))
    data.feo_tot = data.feo;
    data.feo = [];
    return
elseif any(strcmp(data.Properties.VariableNames,'fe2o3')) & ...
    ~any(strcmp(data.Properties.VariableNames,'feo')) & ...
    ~any(strcmp(data.Properties.VariableNames,'feo_tot')) & ...
    ~any(strcmp(data.Properties.VariableNames,'fe2o3_tot'))
    data.feo_tot = data.fe2o3*fe_factor;
    data.fe2o3 = [];
    return
elseif any(strcmp(data.Properties.VariableNames,'fe2o3_tot')) & ...
    ~any(strcmp(data.Properties.VariableNames,'feo')) & ...
    ~any(strcmp(data.Properties.VariableNames,'feo_tot')) & ...
    ~any(strcmp(data.Properties.VariableNames,'fe2o3'))
    data.feo_tot = data.fe2o3_tot*fe_factor;
    data.fe2o3_tot = [];
    return
end

if ~any(strcmp(data.Properties.VariableNames,'feo_tot'))
    data.feo_tot = nan([height(data) 1]);
end
if ~any(strcmp(data.Properties.VariableNames,'fe2o3_tot'))
    data.fe2o3_tot = nan([height(data) 1]);
end


% Convert FeO to FeOT when Fe2O3 = NaN
ind = (~isnan(data.feo) & data.feo > 0) & ...
    isnan(data.feo_tot) & ...
    (isnan(data.fe2o3) | data.fe2o3 <= 0);
data.feo_tot(ind) = data.feo(ind);

% Convert Fe2O3 + Fe0 to FeOT
ind = ~isnan(data.fe2o3) & ...
    ~isnan(data.feo) & ...
    data.feo > 0 & data.fe2o3 > 0;
data.feo_tot(ind) = data.feo(ind) + data.fe2o3(ind)*fe_factor;

% Fe ratio
data.fe2_fe_tot = nan([height(data) 1]);
data.fe2_fe_tot(ind) = data.feo(ind)/molecularwt('FeO') ...
    ./ (data.feo(ind)/molecularwt('FeO') + 2*data.fe2o3(ind)/molecularwt('Fe2O3'));

% Add Fe2O3 to Fe2O3T for Fe2O3T, FeO = NaN
ind = (~isnan(data.fe2o3) & data.fe2o3 > 0) & ...
    isnan(data.fe2o3_tot) & ...
    (isnan(data.feo) | data.feo <= 0);
data.fe2o3_tot(ind) = data.fe2o3(ind);

% Convert Fe2O3T to FeOT for FeOT = NaN
ind = ~isnan(data.fe2o3_tot) & ...
    isnan(data.feo_tot);
data.feo_tot(ind) = data.fe2o3_tot(ind)*fe_factor;

% check for issues with Fe conversion
%d = abs(data.feo_tot - fe_factor*data.fe2o3_tot);
%plot(data.feo_tot,d,'.');
%sum(d > 0.1 & ~isnan(d))
%sum(d > 0.01 & ~isnan(d))
%sum(d > 0.001 & ~isnan(d))
%ind = d > 0.03;
%issues(data,ind)



% remove other iron columns
data.feo = [];
data.fe2o3 = [];
data.fe2o3_tot = [];

return
