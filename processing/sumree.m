function data = sumree(data,scheme)
% SUMREE - Sums rare earth elements
%
%   [ree,lree,mree,hree] = sumree(data,scheme) produces sums of rare-earth
%   elements (REE), as well as subgroups light (LREE), middle (MREE) and
%   heavy (HREE) rare earth elements.  The sums are produced only using
%   values above detections which may have the potential to bias the sums
%   towards higher values when there are large numbers of below detection
%   limit (bdl) samples.
%
%   The scheme has two options:
%       'full'      a sum of all the REE's [La, Ce, Pr, Nd, Sm, Eu, Gd, Tb,
%                   Dy, Ho, Er, Tm, Yb, and Lu] where LREE is defined as
%                   La to Gd and HREE is defined as Tb to Lu
%       'reduced'   a reduced sum to reduce the number of bdl values, in
%                   this case the sum of REE is [La, Ce, Nd, Sm, Eu, Gd,
%                   Er, Yb, and Lu], LREE is [La, Ce, and Nd], HREE is [Er,
%                   Yb, and Lu]
%
%   In both schemes, MREE is [Sm, Eu, and Gd]

% select elements for summation schemes
mree = {'sm_ppm','eu_ppm','gd_ppm'};
switch scheme
    case 'full'
        lree = {'la_ppm','ce_ppm','pr_ppm','nd_ppm','sm_ppm','eu_ppm','gd_ppm'};
        hree = {'tb_ppm','dy_ppm','ho_ppm','er_ppm','tm_ppm','yb_ppm','lu_ppm'};
        ree = [lree, hree];
    case 'reduced'
        lree = {'la_ppm','ce_ppm','nd_ppm'};
        hree = {'er_ppm','yb_ppm','lu_ppm'};
        ree = [lree, mree, hree];
    otherwise
        error('Unknown REE scheme.');
end

% sum REEs
data.ree = sum(data{:,ree},2);      % all REE
data.lree = sum(data{:,lree},2);    % light REE
data.mree = sum(data{:,mree},2);    % middle REE
data.hree = sum(data{:,hree},2);    % heavy REE

% set all REE sums with bdl values as nan
ind = all(data{:,ree} > 0);
data.ree(ind) = nan;
ind = all(data{:,lree} > 0);
data.lree(ind) = nan;
ind = all(data{:,hree} > 0);
data.mree(ind) = nan;
ind = all(data{:,mree} > 0);
data.hee(ind) = nan;

return