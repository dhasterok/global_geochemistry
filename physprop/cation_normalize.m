% Normalizes by moles of cations
function t = cation_normalize(t,varargin)
% CATION_NORMALIZE - Normalizes masses to cations
%
%   t = cation_normalize(t,fact) normalizes the cations in the table t to
%   a value fact, which is optional.  If fact is not supplied, a default
%   value of 100 is used.
%
%   because the table may include a bunch of unused fields,
%   cation_normalize uses global variables CATIONS, OXIDES, NCAT and MW,
%   for the cation field names, oxide field names, number of cations in the
%   oxide, and molecular weight.

global CATIONS OXIDES NCAT MW

% normalization factor
fact = 100;
if nargin == 2
    fact = varargin{1};
end

nr = height(t);
% convert to cations
t{:,CATIONS} = repmat(NCAT,nr,1) .* t{:,lower(OXIDES)} ./ repmat(MW,nr,1);
% now normalize to fact
t{:,CATIONS} = fact * t{:,CATIONS} ./ sum(t{:,CATIONS},2);

return