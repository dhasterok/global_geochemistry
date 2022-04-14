function data = fefix2(data);

%PUT IN CHECKS - TELLING HOW MANY DATA REMOVED VIA QUESTIONABLE METHODS

temp = data(:,{'feo','fe2o3','feo_tot','fe2o3_tot'});

% conversion factor for Fe2O3 to FeO
fe_factor = 2*molecularwt('FeO')/molecularwt('Fe2O3');

%Cases:
temp.fe_fe2_ratio = nan([height(temp) 1]);
temp.feo_tot_c = nan([height(temp) 1]);
temp.fe2o3_tot_c = nan([height(temp) 1]);

%NaN all >100 values
ind = temp.feo > 100;
temp.feo(ind) = NaN;
ind = temp.fe2o3 > 100;
temp.fe2o3(ind) = NaN;
ind = temp.feo_tot > 100;
temp.feo_tot(ind) = NaN;
ind = temp.fe2o3_tot > 100;
temp.fe2o3_tot(ind) = NaN;

%Calculated values - looking at feo and fe2o3 only
%--------------------------------------------------------------------------
%1. feo alone
ind = (~isnan(temp.feo) & isnan(temp.fe2o3));
temp.feo_tot_c(ind) = temp.feo(ind);

%2. fe2o3 alone
ind = (~isnan(temp.fe2o3) & isnan(temp.feo));
temp.fe2o3_tot_c(ind) = temp.fe2o3(ind);

%3. Both feo and fe2o3 and positive, + calculate ratio
ind = (~isnan(temp.feo) & temp.feo >= 0 & ~isnan(temp.fe2o3) & temp.fe2o3 >= 0);
temp.feo_tot_c(ind) = temp.feo(ind) + fe_factor.*temp.fe2o3(ind);
temp.fe_fe2_ratio(ind) = temp.feo(ind)./(temp.feo(ind) + fe_factor.*temp.fe2o3(ind));

%3.1 FeO positive, fe2o3 negative
ind = (~isnan(temp.feo) & temp.feo >= 0 & ~isnan(temp.fe2o3));
temp.feo_tot_c(ind) = temp.feo(ind);

%3.2 Fe2o3 positive, feo negative
ind = (~isnan(temp.fe2o3) & temp.fe2o3 >= 0 & ...
    ~isnan(temp.feo) & temp.feo<0);
temp.fe2o3_tot_c(ind) = temp.fe2o3(ind);

%3.3 Both negative
ind = (~isnan(temp.fe2o3) & ~isnan(temp.feo) & temp.feo < 0 & temp.fe2o3 < 0);
temp.feo_tot_c(ind) = temp.feo(ind);


temp.feo = [];
temp.fe2o3 = [];

%Existing values - corrections/removal
%--------------------------------------------------------------------------

%6. Compare/replace feo_tot and fe2o3_tot if inconsistent with each other
%Both positive but inconsistent with each other

%NOTE: Can actually move these to feo/fe2o3 at some step, as its likely.
ind = (~isnan(temp.feo_tot) & ~isnan(temp.fe2o3_tot) & ...
    temp.feo_tot >= 0 & temp.fe2o3_tot >= 0 & ...
    abs((fe_factor.*temp.fe2o3_tot - temp.feo_tot)./(temp.feo_tot)) > 0.03);
temp.feo_tot(ind)=NaN;
temp.fe2o3_tot(ind)=NaN;
warning_var = size(find(ind==1),1);
fprintf('NaN-ing %d entries: Inconsistent feo & fe2o3 totals\n',warning_var)

%6.2 fe2o3_tot alone
ind = (~isnan(temp.fe2o3_tot) & isnan(temp.feo_tot));
temp.feo_tot(ind) = temp.fe2o3_tot(ind).*fe_factor;

%7. If feo_tot negative and fe2o3_tot positive
ind = (~isnan(temp.fe2o3_tot) & ~isnan(temp.feo_tot) & temp.feo_tot < 0 & temp.fe2o3_tot >= 0);
temp.feo_tot(ind) = temp.fe2o3_tot(ind).*fe_factor;


temp.fe2o3_tot = [];


%Compare calculated and existing: feo_tot_c and feo_tot
%--------------------------------------------------------------------------
%If the difference is greater than 3% - NaN both. OPTIONAL.
%8. Move fe2o3_tot_c if its only one, into feo_tot_c
ind = (~isnan(temp.fe2o3_tot_c) & (isnan(temp.feo_tot_c) | temp.feo_tot_c<0));
temp.feo_tot_c(ind) = temp.fe2o3_tot_c(ind).*fe_factor;
temp.fe2o3_tot_c = [];

%9. Compare consistency between feo_tot and feo_tot_c where exist and
%positive - NaN if inconsistent
ind = (abs((temp.feo_tot_c - temp.feo_tot)./(temp.feo_tot)) > 0.03 & ...
    temp.feo_tot_c >= 0 & temp.feo_tot >= 0);


% Check why its wrong
T = data;
a = temp.feo_tot_c;
T = addvars(T,a,'After','feo_tot');
T(ind,{'feo','fe2o3','fe2o3_tot','feo_tot','a'})




temp.feo_tot_c(ind) = NaN;
temp.feo_tot(ind) = NaN;
warning_var = size(find(ind==1),1);
fprintf('NaN-ing %d entries: Inconsistent feo_tot & feo_tot_c\n',warning_var)

%If feo_tot_c positive but feo_tot negative, replace it
ind = (~isnan(temp.feo_tot_c) & ~isnan(temp.feo_tot) & temp.feo_tot < 0 ...
    & temp.feo_tot_c >= 0);
temp.feo_tot(ind) = temp.feo_tot_c(ind);

%If only feo_tot_c exists, move it
ind = (~isnan(temp.feo_tot_c) & isnan(temp.feo_tot));
temp.feo_tot(ind) = temp.feo_tot_c(ind);

temp.feo_tot_c = [];
data.feo = [];
data.fe2o3 = [];
data.feo_tot = temp.feo_tot;
data.fe2o3_tot = [];
data.fe_fe2_ratio = temp.fe_fe2_ratio;

return