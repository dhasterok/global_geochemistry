function el = fixfieldnames(fields,el)
% FIXFIELDNAMES - fixes field names to match data fields
%
%   el = fixfieldnames(fields,el) fixes the input element or text to match
%   a field name.
if iscell(el)
    for i = 1:length(el)
        el{i} = fixfieldnames(fields,el{i});
    end
else
    % convert spaces to underscores
    ind = strfind(el,' ');
    el(ind) = '_';
    if any(strcmp(fields,el))
        el = el;
    elseif any(strcmp(fields,lower(el)))
        el = lower(el);
    elseif any(strcmpi(fields,[el,'_ppm']))
        el = [lower(el),'_ppm'];
    else
        error(['Could not find field (',el,').']);
    end
end

return