combs = unique(sort(nchoosek(repmat([1 2 3 4],1,4),4),2),'rows');
perms = [];
for i = 1:length(combs)
    perms_comb = uperm(combs(i,:));
    perms = vertcat(perms,perms_comb);
end

perms(perms == 1) = NaN;
perms(perms == 2) = -5;
perms(perms == 3) = 15;
perms(perms == 4) = 127;

T = table(perms(:,1),perms(:,2),perms(:,3),perms(:,4));
T.Properties.VariableNames = {'feo','fe2o3','feo_tot','fe2o3_tot'};
Told = T;

F = fefix(T);
F.Properties.VariableNames = {'feo_calc_tot','fe2_fe_tot'};

C = fefix2(T);
C.Properties.VariableNames = {'feo_calc_tot_2','fe2_fe_tot_2'};

T = [Told F C];
T = sortrows(T);