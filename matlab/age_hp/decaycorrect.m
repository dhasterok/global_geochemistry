function CC = decaycorrect(el,C,age,ref)
% DECAYCORRECT - decay corrects radiogenic isotopes
%
%   CC = decaycorrect(el,C,age,ref) where el is a string with the element
%   to correct, C is the concentration of the element, age is the time in
%   the past in Ma, and ref is the reference for constants.
%
%   At the moment only parent concentrations are corrected.  Parents
%   available are K, Rb, Sm, Th, and U.
%
%   References for constants include:
%      HG17 - in house values (unpublished)
%      R88 - Rybach (1988)
%      TS14 - Turcotte and Schubert (2014)
%      D12 - Dye (2012)
%
%   Note values for Rb and Sm are the same for all references as only HG17
%   addresses these.

% 30 Jul 2018 by D. Hasterok (Univ. Adelaide)

% ensure column vectors
C = C(:);
if length(age) == 1
    age = age*ones(size(C));
else
    age = age(:);
end

% number of samples
ns = length(C);

if length(el > 2)
    el = el(1:end-4);
end

% load isotope constants for correction
const = readtable('radconst.csv','Format','%s%s%f%f%f%f%f');
% limit constants to chosen reference and element corrected
ind = strcmpi(const.reference,ref) & strcmpi(const.element,el);
const = const(ind,:);

% number of radiogenic isotopes
ni = height(const);

% stable isotope concentration
Cs = (1 - sum(const.abundance))*C;

% radiogenic isotope concentration
Cr = repmat(const.abundance',ns,1).*repmat(C,1,ni);
arg = log(2)*repmat(age,1,ni)./repmat(const.half_life',ns,1);

% compute initial concentration at age
CC = Cs + sum(Cr.*exp(arg),2);

return