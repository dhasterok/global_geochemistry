function R = radarprep(data,typefield,types,fields,varargin)
% RADARPLOT - produces a radar/spider plot
%
%   radarprep(DATA,TYPEFIELD,TYPES,FIELDS) where DATA is a table containing
%   both the type field used (TYPEFIELD) to partition the data by TYPES.
%   The spokes of the web are listed in the FIELDS.  Both TYPEFIELD and
%   FIELDS are cell arrays.

vals = zeros([length(types) length(fields)]);
for i = 1:length(types)
    ind = strcmpi(data{:,typefield},types{i});
    vals(i,:) = nanmean(data{ind,fields},1);
    %for j = 1:length(vars)
%         quant_min(i,j) = quantile(data{rock,vars(j)},0.05);
%         quant_min(i,j) = nanmin(data{rock,vars(j)});
%         quant_max(i,j) = quantile(data{rock,vars(j)},0.95);
    %end
end

R.fields = fields;
R.types = types;
R.vals = vals;

return