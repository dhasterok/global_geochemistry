function country_alpha_code = country_code(country_list)
% Takes in a country list - if empty or unknown, ignore it, otherwise try 
% and find it within country_codes.xlsx

country_list = upper(country_list);
country_alpha_code = country_list;
%country_list = cellfun(@lower,country_list,'UniformOutput',0);
[~,~,raw] = xlsread('country_codes.xlsx');
raw = raw(2:end,:);
r = cellfun(@(s)find(strcmpi(s,raw(:,1))),country_list,'uni',0);
r_notempty = ~cellfun(@isempty,r);

if length(find(~r_notempty)) > 0
    not_valid = unique(country_list(find(~r_notempty)));
    for i = 1:size(not_valid,1)
        fprintf('%s\n',not_valid{i})
    end
    fprintf('are not valid countries\n')
    
    
%     for i = 1:length(find(~r_notempty))
%         a = find(~r_notempty);
%         fprintf('%s, ',country_list{a(i)})
%     end
%     fprintf('are not valid countries\n')
end

country_alpha_code(r_notempty) = upper(raw(cell2mat(r),2));

return

