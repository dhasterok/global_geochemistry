function data = cat2ox(data)

% list of cations and their oxide to convert
list = {'Si','SiO2'; 'Ti','TiO2'; 'Al','Al2O3'; ...
    'Cr','Cr2O3'; 'Fe','FeO'; 'Mg','MgO'; 'Ca','CaO'; ...
    'Mn','MnO'; 'Ni','NiO'; 'Ba','BaO'; 'Sr','SrO'; ...
    'Na','Na2O'; 'K','K2O'; 'P','P2O5'};

for i = 1:length(list)
    % field names
    el = [lower(list{i,1}),'_ppm'];
    ox = lower(list{i,2});
    
    if all(~strcmp(el,data.Properties.VariableNames))
        fprintf('  no data found for %s.\n',el);
        continue;
    end
    
    if all(~strcmp(ox,data.Properties.VariableNames))
        data{:,{ox}} = nan([height(data) 1]);
    end
    
    % convert ppm to oxide
    oxide = convert_cat2ox(list(i,:),data{:,{el}});
    
    if strcmp(ox,'feo')
        ox = [ox,'_tot'];
    end
    
    % assign cation converted data where no data exist.
    %a = sum(data{:,ox} > 0);
    ind = oxide > 0 & ~(data{:,{ox}} > 0);
    
    % TEMPORARY CHECK:
%         % Pause
%         fprintf('%s\n',ox)
%         ind2 = ind & oxide > 25;
%         data{ind,{ox}} = oxide(ind);
%         data(ind2,{'filename','sample_name',ox,el})
%         x = input('\n');
    
    
    data{ind,{ox}} = oxide(ind);
    %b = sum(data{:,ox} > 0);
    %[a,b]
end

return

function oxide = convert_cat2ox(field,el,ox)

% determine multiplier
n = 1;
[m,ok] = str2num(field{2}(length(field{1})+1));
if ok
    n = m;
end

f = molecularwt(field{2}) / (1e4 * n * molecularwt(field{1}));

% convert cation to oxide
oxide = f * el;

% remove clearly incorrect values
oxide(oxide > 100) = nan;

return