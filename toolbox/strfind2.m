function ind = strfind2(c,s)
% STRFIND2 - logical strfind
%
%  ind = strfind2(c,s), when all you want is to know whether an item of a
%  cell array (c) contains the string array (s).

list = strfind(c,s);

ind = logical(zeros(size(list)));

for i = 1:length(list)
    ind(i) = isempty(list{i});
end

return