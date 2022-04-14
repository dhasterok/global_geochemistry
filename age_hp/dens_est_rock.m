function data_save = dens_est_rock(data)
%Estimates density for all combinations of major elements + chemical
%indices for a set of rock types
%Rock types to use:
rock_names = {...
    'alkalic basalt','alkalic gabbro';...
    'subalkalic basalt','subalkalic gabbro';...
    'basaltic andesite','gabbroic diorite';...
    'andesite','diorite';...
    'dacite','granodiorite';...
    'rhyolite','granite';...
    'trachybasalt','monzogabbro';...
    'basaltic trachyandesite','monzonite';...
    'trachyandesite','monzonite';...
    'trachydacite','quartz monzonite';...
    'tephrite','foid gabbro';...
    };

major_ele = {...
    'sio2','tio2','al2o3','feo_tot','mgo','cao','k2o','na2o','p2o5','loi',...
    'MALI','ASI','maficity'...
    };

%Fix bad densities
data.density(data.density > 0 & data.density < 4,:) = data.density(data.density > 0 & data.density < 4,:)*1000;

%Remove wrong measurements
data = data(data.sio2 >= 0 & data.sio2 < 100,:);

%Make sure all positive and within range
for i = 1:length(major_ele)
    %If negative set to 0
    data.(major_ele{i})(data.(major_ele{i}) < 0,:) = 0;
    
    %Check if smaller than 100 OR can be NaN
    data = data(data.(major_ele{i}) < 100 | isnan(data.(major_ele{i})),:);
end


density_table = table(rock_names);
density_table.Properties.VariableNames = {'rock_name'};

for i = 1:size(rock_names,1)
    ind = ~isnan(data.density) & (strcmpi(data.rock_type,rock_names{i,1}) | strcmpi(data.rock_type,rock_names{i,2}));
    density_table.ind{i} = ind;
    density_table.num_val(i) = sum(ind);
end

%Density estimation begin
%--------------------------------------------------------------------------
data_save = data;
data_save.density_matt(:,1) = 0;

warning('off','all')

for i = 1:size(rock_names,1)
    
data_save_ind = density_table.ind{i};
data = data_save(data_save_ind,:);
    
y = data.density;
for i = 1:length(major_ele)
    x{:,i} = data.(major_ele{i});
end

n_var = size(x,2);
comb_cell = cell(1,n_var-1);
h = waitbar(0,'Initializing waitbar...');
for i = 1:n_var-1
    %Combinations must include sio2
    c = [ones(size(combnk(2:n_var,i),1),1) combnk(2:n_var,i)];
    b = [];
    for j = 1:size(c,1)
        ind = [];
        for k = 1:size(c,2)
            ind{k} = find(~isnan(x{:,c(j,k)}));
        end
        
        com_ind = ind{1};
        for k = 2:size(c,2)
            com_ind = intersect(com_ind,ind{k});
        end

        ind = com_ind;

        if ~isempty(ind)
            X = [ones(size(x{:,1}))];
            for k = 1:size(c,2)
                X = [X x{:,c(j,k)}];
            end
            b(j,:) = regress(y,X)';
        else
            continue;
        end

        perc = 100*(i/(n_var-1)) + 10*(j/size(c,1));
        waitbar(perc/100,h,sprintf('Progress \n %0.2f%%',perc))
    end
    comb_cell{i} = {c b};
end
close(h)
output = comb_cell;
clearvars  X x y data_dens comb_cell

temp = data(:,{'density','sio2','tio2','al2o3','feo_tot','mgo','cao','k2o','na2o','p2o5','loi',...
    'MALI','ASI','maficity'});
temp.density_matt(:,1) = 0;
h = waitbar(0,'Initializing waitbar...');
for i = 1:size(temp,1)
    sum_nan = sum(~isnan(temp{i,3:end-2}));
        
    if sum_nan == 0
        continue;
    else
        if isnan(temp.sio2(i))
            continue;
        else
            if size(output{sum_nan}{1,1},1)==0
                continue;
            end
            if length(output{sum_nan}{1,2})<2
                continue;
            end
            for j = 1:size(output{sum_nan}{1,1},1)
                if isequal(find(~isnan(temp{i,2:end-2})),output{sum_nan}{1,1}(j,:))
                    %SET DENSITY
                    ind_b = output{sum_nan}{1,1}(j,:);
                    b = output{sum_nan}{1,2}(j,:);
                    temp.density_matt(i,1) = b(1);

                    for k = 2:length(b)
                        temp.density_matt(i,1) = temp.density_matt(i,1) + b(k)*temp{i,ind_b(k-1)+1};
                    end
                    break;
                end
            end
        end
    end
    perc = 100*(i/(size(temp,1)));
    waitbar(perc/100,h,sprintf('Progress \n %0.2f%%',perc))
end
close(h)
temp.density_matt(temp.density_matt==0) = NaN;

data_save.density_matt(data_save_ind,:) = temp.density_matt;

end
data_save.density_matt(data_save.density_matt==0) = NaN;

data_save = densest(data_save);

warning('on','all')
return