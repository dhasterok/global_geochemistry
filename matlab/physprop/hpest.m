function data = hpest(data,varargin)
% HPEST - estimates heat production from a data table
%
% data = hpest(data) computes the heat production for a table adding
% a column for heat_production_mass (uW kg^-1) and heat_production
% (uW m^-3)


% make sure all potassium is expressed as k2o
k2o = getK2O(data);
th = data.th_ppm;
u = data.u_ppm;

% deal with censored data (values below detection limit)
% default below detection limit values [K2O(wt%), Th(ppm), U(ppm)]
bdl_values = [-0.01 -1 -1];
if nargin == 2
    bdl_values = varargin{1};
end
fprintf('Using the BDL(0) values:\n');
fprintf('  K2O: %f\n',bdl_values(1));
fprintf('  Th: %f\n',bdl_values(2));
fprintf('  U: %f\n',bdl_values(3));

k2o(k2o == 0) = bdl_values(1);
th(th == 0) = bdl_values(2);
u(u == 0) = bdl_values(3);

% sign for indicating one element is below detection limit.
s = nan(size(k2o));
% uncensored data
s(k2o > 0 & th > 0 & u > 0) = 1;
% censored data
s(k2o <= 0 | th <= 0 | u <= 0) = -1;

% heat production by mass
data.heat_production_mass = s.*computehp(abs(k2o),abs(u),abs(th));

if ~any(strcmp('density_model',fieldnames(data)))
    data = vpest(data);
end

data.heat_production = data.heat_production_mass.*data.density_model;

% when density_model cannot be determined, use density_sio2 to estimate
% heat_production
%ind = isnan(data.density_model) & ~isnan(data.heat_production_mass) & data.sio2 >= 40;
%data.heat_production(ind) = data.heat_production_mass(ind).*data.density_sio2(ind);

%rho = nan([length(ind) 1]);
%if isfield(data,'COMPOSITION')
%    for i = 1:length(ind)
%        if strcmp(data.COMPOSITION(ind(i)),'FELSIC')
%            rho(i) = 2700;
%        elseif strcmp(data.COMPOSITION(ind(i)),'INTERMEDIATE')
%            rho(i) = 2850;
%        elseif strcmp(data.COMPOSITION(ind(i)),'MAFIC')
%            rho(i) = 3000;
%        elseif strcmp(data.COMPOSITION(ind(i)),'ULTRAMAFIC')
%            rho(i) = 3300;
%        end
%    end
%else
%    rho = NaN;
%end

%data.heat_production(ind) = rho .* data.heat_production_mass(ind);

return
