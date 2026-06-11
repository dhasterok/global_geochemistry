function spreadsheet_names = database_table(data,numrows)
%Counts number of values in data for each spreadsheet_id and orders them
%Top 10 are listed, and rest as misc. sources.
%postgres_earthchem files are combined
%postgres_geoatlas files are combined
%postgres_of_ files are combined

%Create a cell array containing all filenames and number of values
spreadsheet_names = {};
spreadsheet_names(:,1) = unique(data.filename);
for i = 1:length(spreadsheet_names(:,1))
    spreadsheet_names{i,2} = length(find(strcmpi(data.filename,spreadsheet_names{i,1})));
    spreadsheet_names{i,3} = (find(strcmpi(data.filename,spreadsheet_names{i,1})));
end

%Combine earthchem, geoatlas, and canada files, and remove individual entry
ind_earth = strfind(spreadsheet_names(:,1), 'postgres_earthchem');
ind_earth = find(not(cellfun('isempty', ind_earth)));
ind_atlas = strfind(spreadsheet_names(:,1), 'postgres_geoatlas');
ind_atlas = find(not(cellfun('isempty', ind_atlas)));
ind_canada = strfind(spreadsheet_names(:,1), 'postgres_of_');
ind_canada = find(not(cellfun('isempty', ind_canada)));

ind_OZCHEM = strfind(spreadsheet_names(:,1), 'postgres_geochemistry_AUS');
ind_OZCHEM = find(not(cellfun('isempty', ind_OZCHEM)));
%spreadsheet_names{ind_OZCHEM,1} = 'OZCHEM';

spreadsheet_names{end+1,1} = 'EarthChem';
spreadsheet_names{end,2} = sum(cell2mat(spreadsheet_names(ind_earth,2)));

spreadsheet_names{end,3} = [];
for i = 1:length(ind_earth)
    spreadsheet_names{end,3} = union(spreadsheet_names{end,3},spreadsheet_names{ind_earth(i),3});
end

spreadsheet_names{end+1,1} = 'Geoscience Atlas';
spreadsheet_names{end,2} = sum(cell2mat(spreadsheet_names(ind_atlas,2)));
spreadsheet_names{end,3} = [];
for i = 1:length(ind_atlas)
    spreadsheet_names{end,3} = union(spreadsheet_names{end,3},spreadsheet_names{ind_atlas(i),3});
end

spreadsheet_names{end+1,1} = 'Canadian Database of Geochemical Surveys';
spreadsheet_names{end,2} = sum(cell2mat(spreadsheet_names(ind_canada,2)));
for i = 1:length(ind_canada)
    spreadsheet_names{end,3} = union(spreadsheet_names{end,3},spreadsheet_names{ind_canada(i),3});
end

spreadsheet_names(vertcat(ind_earth,ind_atlas,ind_canada),:) = [];




[trash idx] = sort([spreadsheet_names{:,2}], 'descend');
spreadsheet_names = spreadsheet_names(idx,:);

%Specify number of 'named' databases in the table
n = numrows;
sum_ind = (n+1):1:length(spreadsheet_names(:,1));
sum_val = sum(cell2mat(spreadsheet_names(sum_ind,2)));
short_table = vertcat(spreadsheet_names(1:n,1:2),{'Misc.',sum_val});


% print table of number of rocks in each field
fprintf(' # Data                  Database\n');
fprintf('----------   -------------------------------\n');
for i = 1:size(short_table,1)
    fprintf('%9i    %s\n',short_table{i,2},short_table{i,1});
end





%Create a cell array containing all filenames and number of values
spreadsheet_names2 = {};
spreadsheet_names2(:,1) = unique(data.filename);

%Select ages greater than 200Ma (the first bin)
for i = 1:length(spreadsheet_names2(:,1))
    spreadsheet_names2{i,2} = length(find(strcmpi(data.filename(data.avg_age>=200,:),spreadsheet_names2{i,1})));
end

%Combine earthchem, geoatlas, and canada files, and remove individual entry
ind_earth = strfind(spreadsheet_names2(:,1), 'postgres_earthchem');
ind_earth = find(not(cellfun('isempty', ind_earth)));
ind_atlas = strfind(spreadsheet_names2(:,1), 'postgres_geoatlas');
ind_atlas = find(not(cellfun('isempty', ind_atlas)));
ind_canada = strfind(spreadsheet_names2(:,1), 'postgres_of_');
ind_canada = find(not(cellfun('isempty', ind_canada)));

ind_OZCHEM = strfind(spreadsheet_names2(:,1), 'postgres_geochemistry_AUS');
ind_OZCHEM = find(not(cellfun('isempty', ind_OZCHEM)));
%spreadsheet_names2{ind_OZCHEM,1} = 'OZCHEM';

spreadsheet_names2{end+1,1} = 'EarthChem';
spreadsheet_names2{end,2} = sum(cell2mat(spreadsheet_names2(ind_earth,2)));
spreadsheet_names2{end,3} = [];
for i = 1:length(ind_earth)
    spreadsheet_names2{end,3} = union(spreadsheet_names2{end,3},spreadsheet_names2{ind_earth(i),3});
end
spreadsheet_names2{end+1,1} = 'Geoscience Atlas';
spreadsheet_names2{end,2} = sum(cell2mat(spreadsheet_names2(ind_atlas,2)));
spreadsheet_names2{end,3} = [];
for i = 1:length(ind_atlas)
    spreadsheet_names2{end,3} = union(spreadsheet_names2{end,3},spreadsheet_names2{ind_atlas(i),3});
end
spreadsheet_names2{end+1,1} = 'Canadian Database of Geochemical Surveys';
spreadsheet_names2{end,2} = sum(cell2mat(spreadsheet_names2(ind_canada,2)));
for i = 1:length(ind_canada)
    spreadsheet_names2{end,3} = union(spreadsheet_names2{end,3},spreadsheet_names2{ind_canada(i),3});
end
spreadsheet_names2(vertcat(ind_earth,ind_atlas,ind_canada),:) = [];

[trash idx] = sort([spreadsheet_names2{:,2}], 'descend');
spreadsheet_names2 = spreadsheet_names2(idx,:);

%Specify number of 'named' databases in the table
n = numrows;
sum_ind = (n+1):1:length(spreadsheet_names2(:,1));
sum_val = sum(cell2mat(spreadsheet_names2(sum_ind,2)));
short_table = vertcat(spreadsheet_names2(1:n,1:2),{'Misc.',sum_val});


% print table of number of rocks in each field
fprintf('# Data>200              Database\n');
fprintf('----------   -------------------------------\n');
for i = 1:size(short_table,1)
    fprintf('%9i    %s\n',short_table{i,2},short_table{i,1});
end


return