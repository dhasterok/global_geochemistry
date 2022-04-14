function [T,ind] = avgchem(data,varargin)
% AVGCHEM - Computes average chemistry and statistics for a given rock
% type.
%
%    T = avgchem(table_of_chem,elements,rtype) where rtype is the rock
%    types as a cell array, and elements are the elements or oxides to
%    average.
%
%    Example:
%       T = avgchem(data, ...
%           {'SiO2','TiO2','Al2O3','FeO_tot', ...
%           'MnO','MgO','CaO','Na2O','K2O','P2O5'}, ...
%            {'dolomite','dolomitic','dolostone'});

% determine indicies with rock types in rtype
if nargin > 1 & mod(nargin,2) ~= 1
    error('Input option not defined.');
end

% default chemical output
elements = {'SiO2','TiO2','Al2O3','FeO','MnO', ...
    'MgO','CaO','Na2O','K2O','P2O5', ...
    'Sc','V','Cr','Co','Ni','Cu','Zn','As', ...
    'Rb','Sr','Y','Zr','Nb','Ba','La','Ce', ...
    'Pr','Nd','Sm','Eu','Gd','Tb','Dy','Ho', ...
    'Er','Tm','Yb','Lu','Hf','Ta','W','Pb', ...
    'Th','U'};

opt = 1;
rtype = {};
while opt < nargin
    switch lower(varargin{opt})
    case 'rocktype'
        rtype = varargin{opt+1};
    case 'elements'
        elements = varargin{opt+1};
    end
    opt = opt + 2;
end

if ~isempty(rtype)
    ind = logical(zeros([height(data) 1]));
    for k = 1:length(rtype)
        tmp = strfind(lower(data.rock_name),rtype{k});
        for i = 1:length(tmp)
            if ~isempty(tmp{i})
                ind(i) = 1;
            end
        end
    end
else
    ind = logical(ones([height(data) 1]));
end

% compute stats on given elements
for i = 1:length(elements)
    if any(strcmpi(lower(elements{i}),data.Properties.VariableNames))
        if strcmp(lower(elements{i}),'feo')
            [t(i,1),t(i,2),t(i,3),t(i,4:8),t(i,9),t(i,10)] = ...
                statsnonan(data{ind,[lower(elements{i}),'_tot']});
        else
            [t(i,1),t(i,2),t(i,3),t(i,4:8),t(i,9),t(i,10)] = ...
                statsnonan(data{ind,lower(elements{i})});
        end
    elseif any(strcmpi([lower(elements{i}),'_ppm'],data.Properties.VariableNames))
        [t(i,1),t(i,2),t(i,3),t(i,4:8),t(i,9),t(i,10)] = ...
            statsnonan(data{ind,[lower(elements{i}),'_ppm']});
    else
        [lower(elements{i}),'_tot']
    end
end

% convert array to table
T = array2table(t,'VariableNames',{'N','mu','sigma','Q025','Q250','Q500','Q750','Q975','min','max'});

% add element names to table
T.elements = elements(:);

return

function [n,v_avg,v_sd,v_q,v_min,v_max] = statsnonan(v)

q = [0.025 0.25 0.5 0.75 0.975];
v = v(v > 0);

n = length(v);
if n == 0
    v_avg = NaN;
    v_sd = NaN;
    v_min = NaN;
    v_max = NaN;
    v_q = nan(size(q));
    return
end
v_avg = mean(v);
v_sd = std(v);
v_q = quantile(v,q);
v_min = min(v);
v_max = max(v);

return
