function [maintable] = merge_table(maintable,subtable)
% Merge tables based on ID - automatically found
% Merges split subtables back into one table - must share a KEY column e.g.
% sample_id

% Find ID column
id_column = subtable.Properties.VariableNames(ismember(...
    subtable.Properties.VariableNames,maintable.Properties.VariableNames));

if length(id_column) > 1
    error('Duplicate columns/more than one ID column?')
    return
end

[~,subtable_ind] = link_column_indices(maintable(:,id_column),subtable(:,id_column));
columns_add = subtable.Properties.VariableNames(~ismember(...
    subtable.Properties.VariableNames,maintable.Properties.VariableNames));
maintable = [maintable subtable(subtable_ind,columns_add)];




return