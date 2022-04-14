function [age_min,avg_age,age_max] = age_adjust(age_min,age,age_max,age_sd,age_var)
% avg_age = age_adjust(age_min,age,age_max,age_var)
% Returns: avg_age - reordered/corrected age if necessary, and average of
% min_age and max_age where age did not previously exist.
% age_min:

%avg_age = NaN(size(age));
avg_age = age;

initial_ages = length(find(~isnan(avg_age)));


fprintf('\n----------------------\n')
fprintf('age_adjust corrections\n')
fprintf('----------------------\n\n')

%--------------------------------------------------------------------------
% 1. Removing negatives
%--------------------------------------------------------------------------

% Negative ages: Currently all are less than 10 Ma. For now - remove.
ind = (avg_age < 0);
avg_age(ind) = NaN;
fprintf('age < 0 removed: %d\n',length(find(ind)))

ind = (age_min < 0);
age_min(ind) = NaN;
fprintf('age_min < 0 removed: %d\n',length(find(ind)))

ind = (age_max < 0);
age_max(ind) = NaN;
fprintf('age_max < 0 removed: %d\n',length(find(ind)))

% Ages that are greater than age of Earth - assume in years rather than Ma
ind = (avg_age>4600);
avg_age(ind) = avg_age(ind)./1000000;
fprintf('age > 4600 Ma corrected to years: %d\n',length(find(ind)))

ind = (age_min>4600);
age_min(ind) = age_min(ind)./1000000;
fprintf('age_min > 4600 Ma corrected to years: %d\n',length(find(ind)))

ind = (age_max>4600);
age_max(ind) = age_max(ind)./1000000;
fprintf('age_max > 4600 Ma corrected to years: %d\n',length(find(ind)))

%--------------------------------------------------------------------------
% 2. Add age_min, age_max where age and age_sd exists
%--------------------------------------------------------------------------

% Number of age_sd values
ind = ~isnan(age_sd);

% Number of age and age_sd values
ind2 = ~isnan(age_sd) & ~isnan(avg_age);
fprintf('Age_sd and age not aligned correctly and removed: %d\n',length(find(ind))-length(find(ind2)))

% Adding age_min and age_max to ages with sd only
age_min(ind2) = avg_age(ind2) - age_sd(ind2);
age_max(ind2) = avg_age(ind2) + age_sd(ind2);
fprintf('Added age_min and age_max from age_sd: %d\n',length(find(ind2)))

%--------------------------------------------------------------------------
% 3. Reordering incorrect min/max/avg ages
%--------------------------------------------------------------------------

% Multiple cases for missing age_min, age, or age_max combinations:
% Case 1: Only age_max - could be correct. Leave as is.


% Case 2: Only age_min - same as above


% Case 3: Only age - already copied in


% Case 4: age_min and age_max
% First check/rotate age_min and age_max if wrong
ind = (isnan(avg_age) & (age_min > age_max));
temp = age_min(ind);
age_min(ind) = age_max(ind);
age_max(ind) = temp;
clear temp

% Set age to midpoint between them
ind = (isnan(avg_age) & ~isnan(age_min) & ~isnan(age_max));
avg_age(ind) = (age_max(ind) + age_min(ind))./2;
fprintf('Added avg_age using age_min and age_max average: %d\n',length(find(ind)))


% Case 5: age_min and age
% First check/rotate age_min and age if wrong
ind = (isnan(age_max) & (age_min > avg_age));
temp = age_min(ind);
age_min(ind) = avg_age(ind);
avg_age(ind) = temp;
clear temp

% Case 6: age_max and age
ind = (isnan(age_min) & (avg_age > age_max));
temp = age_max(ind);
age_max(ind) = avg_age(ind);
avg_age(ind) = temp;
clear temp

% Case 7: All 3 values
ind = (~isnan(age_min) & avg_age > age_max);
temp = avg_age(ind);
avg_age(ind) = age_max(ind);
age_max(ind) = temp;
clear temp

ind = (~isnan(avg_age) & age_min > age_max);
temp = age_min(ind);
age_min(ind) = age_max(ind);
age_max(ind) = temp;
clear temp

ind = (~isnan(age_max) & age_min > avg_age);
temp = avg_age(ind);
avg_age(ind) = age_min(ind);
age_min(ind) = temp;
clear temp

%--------------------------------------------------------------------------
% 3. Removing points which dont lie within designated age variance
%--------------------------------------------------------------------------

ind = (age_max - age_min) > age_var;
avg_age(ind) = NaN;

fprintf('Removed avg_age values outside %d Ma tolerance: %d\n',...
    age_var,length(find(ind)))

ind = (isnan(age_max) | isnan(age_min)) & ~isnan(avg_age);
fprintf('Values with unknown/ambiguous range (EDIT: not removed): %d\n',length(find(ind)))
%avg_age(ind) = NaN;

final_ages = length(find(~isnan(avg_age)));
fprintf('No. ages initial: %d\nNo. ages final: %d\n',initial_ages,final_ages)

return