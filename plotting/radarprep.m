function R = radarprep(data,varargin)
% radarprep - prepares a data table for a radar plot.
%
%   radarprep(DATA,TYPEFIELD,TYPES,FIELDS) where DATA is a table containing
%   both the type field used (TYPEFIELD) to partition the data by TYPES.
%   The spokes of the web are listed in the FIELDS.  Both TYPEFIELD and
%   FIELDS are cell arrays.
%
%   An optional argument of quantile levels can be added.
%   radarprep(DATA,TYPEFIELD,TYPES,FIELDS,QLEVELS)

p = inputParser;

addRequired(p,'Data',@istable);
addParameter(p,'Fields',{},@iscell);
addParameter(p,'GroupField','',@ischar);
addParameter(p,'Groups',[]);
addParameter(p,'Quantiles',[],@isnumeric);

parse(p,data,varargin{:});

fields = p.Results.Fields;
groupfield = p.Results.GroupField;
groups = p.Results.Groups;
quanta = p.Results.Quantiles;

if ~isempty(quanta)
    if isempty(groupfield)
        vals = zeros([1 length(fields) length(quanta)]);
    else
        vals = zeros([length(groups) length(fields) length(quanta)]);
    end
else
    vals = zeros([length(groups) length(fields)]);
end

if isempty(groupfield)
    if isempty(quanta)
        vals = nanmean(data{:,fields},1);
    else
        for j = 1:length(quanta)
             vals(1,:,j) = quantile(data{:,fields},quanta(j));
        end
    end
else
    for i = 1:length(groups)
        if iscell(groups)
            ind = strcmpi(data{:,groupfield},groups{i});
        else
            ind = data{:,groupfield} == groups(i);
        end
        if isempty(quanta)
            vals(i,:) = nanmean(data{ind,fields},1);
        else
            for j = 1:length(quanta)
                 vals(i,:,j) = quantile(data{ind,fields},quanta(j));
            end
        end
    end
end

R.fields = fields;
R.types = groups;
R.vals = vals;
R.quanta = quanta;

return