function data = cat_mol_norm(data,oxides)

total = zeros([height(data) 1]);
for i = 1:length(oxides)
    if strcmp(oxides{i},'FeO_tot')
        [str,res] = strtok('FeO','O');
    else
        [str,res] = strtok(oxides{i},'O');
    end

    if length(res) == 1
       res = 1;
    else
        res = str2num(res(end));
    end

    [n,ok] = str2num(str(end));
    if ~ok | imag(n)
        n = 1;
    else
        str = str(1:end-1);
    end

    el{i,1} = str;
    el{i,2} = n;
    el{i,3} = res;
   % {oxides{i}, el{i,1}, el{i,2}}

    if strcmp(oxides{i},'FeO_tot')
        data{:,el{i,1}} = data{:,lower(oxides{i})} * el{i,2} / ...
            molecularwt('FeO');
    else
        data{:,el{i,1}} = data{:,lower(oxides{i})} * el{i,2} / ...
            molecularwt(oxides{i});
    end
    total = total + data{:,el{i,1}};
end

for i = 1:length(oxides)
    data{:,el{i,1}} = data{:,el{i,1}}./total;
    ind = data{:,el{i,1}} <= 0;
    data{ind,el{i,1}} = NaN;
end

return
