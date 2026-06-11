function [output_table_ind,input_table_ind] = link_column_indices(output_table_id,input_table_id)
% Inputs: Shared ID for both tables
% Output: An indices vector for the rows in the input_table that line with
% out output table
[~,idx] = ismember(output_table_id,input_table_id);
output_table_ind = ~(idx==0);
idx(idx == 0) = [];
input_table_ind = idx;

return