function [data,subtable] = partition_db(data,column_list,id_name)
% partition_db splits database into relational database structure subtables
% Inputs:
%   data - table variable
%   column_list - column list to partition off from data and send to new
%                 subtable
%   id_name - KEY column name to utilise to relate the two tables e.g.
%             sample_id
%
% Outputs: data - table variable with columns removed. KEY column added to
%                 correlate to new subtable
%          subtable - new subtable with shared key column ID
%

% Replace all NaN with Inf - this is distinct and matchable for unique and
% ismember
temp = data(:,column_list);
missing_ind = ismissing(temp,NaN);
temp = table2cell(temp);
temp(missing_ind) = {Inf};
data(:,column_list) = temp;

% Create the unique subtable and add an id column
subtable = unique(data(:,column_list));
subtable = [array2table([1:size(subtable,1)]','VariableNames',{id_name}) subtable];

% Find the ID that match in the original data table and add the reference
[~,idata] = ismember(data(:,column_list),subtable(:,column_list));
data.(id_name) = nan([height(data) 1]);
data{:,id_name} = idata;
data(:,column_list) = [];


% Replace all Inf with NaN again
temp = subtable(:,column_list);
missing_ind = ismissing(temp,Inf);
temp = table2cell(temp);
temp(missing_ind) = {NaN};
subtable(:,column_list) = temp;


return