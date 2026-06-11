function duplicate_check(data)
% Crude duplicate check - based on sample name only

[unique_samples, ~, v] = unique(data.sample_name);
occurrence_samples = accumarray(v,1);

% Only check for 2 dupes - 3+ is too complicated
ind = find(occurrence_samples == 2);

count = 0;

h = waitbar(0,'Please wait...');

tol_eq = @(x,y,tol) abs(x-y)<tol;


unique_files = unique(data.filename);
count_overlap = zeros(length(unique_files));


for i = 1:length(ind)
    if isempty(unique_samples{ind(i)})
        continue
    else
        ind2 = find(strcmpi(unique_samples{ind(i)},data.sample_name));
        count = count + size(ind2,1) - 1;
        if tol_eq(data.sio2(ind2(1),:),data.sio2(ind2(2),:),1)
            idx_sp1 = find(strcmpi(unique_files,data.filename(ind2(1),:)));
            idx_sp2 = find(strcmpi(unique_files,data.filename(ind2(2),:)));
            count_overlap(idx_sp1,idx_sp2) = count_overlap(idx_sp1,idx_sp2) + 1;
            count_overlap(idx_sp2,idx_sp1) = count_overlap(idx_sp2,idx_sp1) + 1;
        end
    end
    if mod(i,100)
        waitbar(i/length(ind),h,'% done')
    end
end

num_to_find = 200;

fprintf('\n\n')
fprintf('Top %d overlaps:\n',num_to_find)
total_dupes = 0;

count_overlap = tril(count_overlap,-1);
[max_vals,max_idx] = maxk(count_overlap(:),num_to_find);

index_list = [];
for i = 1:num_to_find
    [index_list(i,1),index_list(i,2)] = ind2sub(size(count_overlap),max_idx(i));
end

% for i = 1:num_to_find
%     fprintf('%s and %s share: %d values\n',unique_files{index_list(i,1)},unique_files{index_list(i,2)},max_vals(i))
%     total_dupes = total_dupes + max_vals(i);
% end

i = 1;
while max_vals(i) > 0
    fprintf('%s and %s share: %d values\n',unique_files{index_list(i,1)},unique_files{index_list(i,2)},max_vals(i))
    total_dupes = total_dupes + max_vals(i);
    i = i+1;
end



fprintf('Total dupes: %d in %d samples\n',total_dupes,size(data,1))




return