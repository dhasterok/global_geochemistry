function avg_age = age_correction(data,age_var)
%[avg_age] = age_correction(data,age_var)
%Data is a table format containing age_min, age and age_max columns
%Variance is the maximum age difference between min and max you want to
%allow: e.g. specifying 200 means all data with >200Ma between min and max
%is NaN'd
%Returns: avg_age - reordered/corrected age if necessary, and average of
%min_age and max_age where age did not previously exist.
%sum((strcmpi('US',data.country) | strcmpi('United States',data.country)) & data.heat_production > 0 & data.age > 0)
%% Removing negatives/reordering incorrect min/max/avg ages
    %All negative ages are less than 10Ma (as of 15/05/2017), so set them positive
    data.age = abs(data.age);
    data.age_min = abs(data.age_min);
    data.age_max = abs(data.age_max);
    
    %If ages are greater than age of earth -> assume its in years rather
    %than millions of years. Theres probably more like this - how to check? e.g.
    %1000 year old volcanics may be listed as 1000.
    ind = find(data.age>4600);
    data.age(ind) = data.age(ind)./1000000;
    if isempty(ind)
        ind = [];
    end
    fprintf('age > 4600 corrected: %d\n',length(ind))
    
    ind = find(data.age_min>4600);
    data.age_min(ind) = data.age_min(ind)./1000000;
    if isempty(ind)
        ind = [];
    end
    fprintf('age_min > 4600 corrected: %d\n',length(ind))
    
    ind = find(data.age_max>4600);
    data.age_max(ind) = data.age_max(ind)./1000000;
    if isempty(ind)
        ind = [];
    end
    fprintf('age_max > 4600 corrected: %d\n',length(ind))
    
    %Check for NaN age but age_min, age_max exist
    %Case 1: Only age_max
        %Could be correct - must be younger than X. Leave age blank then as
        %unknown
    %Case 2: Only age_min
        %Same as above - leave age blank
    %Case 3: Only age
        data.avg_age = data.age;
    
    %Case 6: min_age, max_age
        %First: Swap min_age and max age if needed
        %Age has been saved to avg_age already, just check that min_age < age_max
        ind = find(data.age_max < data.age_min);
        %Swap min_age and age
        temp = data.age_max(ind);
        data.age_max(ind) = data.age_min(ind);
        data.age_min(ind) = temp;
        if isempty(ind)
            ind = [];
        end
        fprintf('min_age > max_age corrected: %d\n',length(ind))
    
        %Second: Average them and place in avg_age - only if age doesnt
        %exist already
        ind = find(isnan(data.avg_age) & ~isnan(data.age_min) & ~isnan(data.age_max));
        data.age(ind) = (data.age_max(ind) + data.age_min(ind))./2;
        data.avg_age(ind) = data.age(ind);
    
    %Case 4: min_age and age
        %Age has been saved to avg_age already, just check that min_age < age
        ind = find(~isnan(data.age) & ~isnan(data.age_min) & data.age_min > data.age);
        %Swap min_age and age
        temp = data.age_min(ind);
        data.age_min(ind) = data.age(ind);
        data.age(ind) = temp;
        data.avg_age(ind) = data.age(ind);
        if isempty(ind)
            ind = [];
        end
        fprintf('age < min_age corrected: %d\n',length(ind))
    
    %Case 5: age and max_age
        %Age has been saved to avg_age already, just check that age < age_max
        ind = find(~isnan(data.age) & ~isnan(data.age_max) & data.age_max < data.age);
        %Swap max_age and age
        temp = data.age_max(ind);
        data.age_max(ind) = data.age(ind);
        data.age(ind) = temp;
        data.avg_age(ind) = data.age(ind);
        if isempty(ind)
            ind = [];
        end
        fprintf('age > max_age corrected: %d\n',length(ind))
        
        %% Setting age variance between min/max
        ind = find(~isnan(data.avg_age) & (data.age_max - data.age_min) > age_var);
        ind2 = find((~strcmpi(data.time_period,'') & strcmpi(data.time_period_min,'') & strcmpi(data.time_period_max,'')));
        ind4 = find((strcmpi(data.time_period,'') & strcmpi(data.time_period_min,'') & ~strcmpi(data.time_period_max,'')));
        ind3 = find((strcmpi(data.time_period,'') & ~strcmpi(data.time_period_min,'') & strcmpi(data.time_period_max,'')));
        ind5 = find((~strcmpi(data.time_period,'') & strcmpi(data.time_period_min,'') & strcmpi(data.time_period_max,''))...
            | (strcmpi(data.time_period,'') & strcmpi(data.time_period_min,'') & ~strcmpi(data.time_period_max,''))...
            | (strcmpi(data.time_period,'') & ~strcmpi(data.time_period_min,'') & strcmpi(data.time_period_max,'')));
        ind6 = find(~isnan(data.avg_age) | (data.age_max - data.age_min) < age_var);
        
        data.avg_age(ind) = NaN;
        
        if isempty(ind)
            ind = [];
        end
        fprintf('%d ages with range greater than %d Mya removed\n',length(ind),age_var)
        sum([data.avg_age >= 0 data.age >= 0])    
        fprintf('%d total data with ages\n',length(ind6))
        
        fprintf('%d can potentially be added using time period\n',length(ind5))
    
        if iscell(data.time_period_min)
            fprintf('time_period_min:\n')
            a=unique(data.time_period_min(ind3,:),'stable');
            b=cellfun(@(x) sum(ismember(data.time_period_min(ind3,:),x)),a,'un',0);
            [a b]
        end

        if iscell(data.time_period)
            fprintf('time_period:\n')
            a=unique(data.time_period(ind2,:),'stable');
            b=cellfun(@(x) sum(ismember(data.time_period(ind2,:),x)),a,'un',0);
            [a b]
        end

        if iscell(data.time_period_max)
            fprintf('time_period_max:\n')
            a=unique(data.time_period_max(ind4,:),'stable');
            b=cellfun(@(x) sum(ismember(data.time_period_max(ind4,:),x)),a,'un',0);
            [a b]
        end
        
        avg_age = data.avg_age;
        
     %sum((strcmpi('US',data.country) | strcmpi('United States',data.country)) & data.heat_production > 0 & data.avg_age > 0)
        %THIS SECTION: Before removing those greater than tolerance:
        %Give them time periods. Era?

return