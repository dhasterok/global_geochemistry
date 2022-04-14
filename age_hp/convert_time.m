function data = convert_time(data)

[nr,~] = size(data.age); 

% adjust from geologic time scale to age
tmp = convert_time_vector(data.age_min);
tmp_min = tmp(:,1);
tmp_avg = convert_time_vector(data.age);
tmp = convert_time_vector(data.age_max);
if numel(tmp) == nr
    tmp_max = tmp;
else
    tmp_max = tmp(:,3);
end

tmp = nan([nr 3]);

% set minimum bound
tmp(:,1) = tmp_min;
if numel(tmp_avg) ~= nr
    ind = (isnan(tmp_min) & ~isnan(tmp_avg(:,1)));
    tmp(ind,1) = tmp_avg(ind,1);
end

% set maximum bound
tmp(:,3) = tmp_max;
if numel(tmp_avg) ~= nr
    ind = (isnan(tmp_max) & ~isnan(tmp_avg(:,3)));
    tmp(ind,3) = tmp_avg(ind,3);
end

% set average
if numel(tmp_avg) ~= nr
    tmp(:,2) = tmp_avg(:,2);
else
    tmp(:,2) = tmp_avg;
end
ind = (isnan(tmp(:,2)) & ~isnan(tmp(:,1)) & ~isnan(tmp(:,3)));
tmp(ind,2) = 0.5*(tmp(ind,1) + tmp(ind,3));
data.age_merge = tmp;

return

function age = convert_time_vector(vect);

if isnumeric(vect)
    age = vect;
    return;
end

age = nan([length(vect),3]);
for i = 1:length(vect)
    if isempty(vect{i})
        continue;
    end

    if ischar(vect{i})
        [val,ok] = str2num(vect{i});
        if ~ok
            age(i,:) = age2ts(vect{i});
        else
            if isempty(val)
                continue;
            end
            age(i,2) = val;
        end
    else
        age(i,2) = vect{i};
    end
end

return
