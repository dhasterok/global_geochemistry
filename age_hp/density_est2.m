function output = density_est2(data)

rock_names = {...
    'picrobasalt','peridotgabbro';...
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
    'trachyte','syenite';...
    'tephrite','foid gabbro';...
    'ultramafic foidite','ultramafic foidolite';...
    };


major_ele = {...
    'sio2','tio2','al2o3','feo_tot','mgo','cao','k2o','na2o','p2o5','loi',...
    'MALI','ASI','maficity'...
    };

%Fix bad densities
data.density(data.density > 0 & data.density < 4,:) = data.density(data.density > 0 & data.density < 4,:)*1000;

%Remove wrong measurements
data = data(data.sio2 >= 0 & data.sio2 < 100,:);

%Make sure all positive and within range if exist
for i = 1:length(major_ele)
    %If negative set to 0
    data.(major_ele{i})(data.(major_ele{i}) < 0,:) = 0;
    
    %Check if smaller than 100 OR can be NaN
    data = data(data.(major_ele{i}) < 100 | isnan(data.(major_ele{i})),:);
end

%Save a copy of adjustments data
saved_data = data;


density_table = table(rock_names);
density_table.Properties.VariableNames = {'rock_name'};

for i = 1:size(rock_names,1)
    ind = ~isnan(data.density) & (strcmpi(data.rock_type,rock_names{i,1}) | strcmpi(data.rock_type,rock_names{i,2}));
    density_table.ind{i} = ind;
    density_table.num_val(i) = sum(ind);
    ind_nodens = (strcmpi(data.rock_type,rock_names{i,1}) | strcmpi(data.rock_type,rock_names{i,2}));
    density_table.ind_nodens{i} = ind_nodens;
    density_table.num_val_nodens(i) = sum(ind_nodens);
end

density_table

%Density estimation begin
%--------------------------------------------------------------------------
data_save = data;
data_save.density_matt(:,1) = 0;

warning('off','all')

for g = 1:size(rock_names,1)

fprintf('%s ...\n',density_table.rock_name{g})
data_save_ind = density_table.ind{g};
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
            if size(X,1) >= 10
                b(j,:) = regress(y,X)';
            else
                b(j,:) = zeros(size(X,2),1);
                sprintf('At least one model is under done \n')
                continue;
            end
        else
            continue;
        end

    end
    perc = 100*(i/(n_var-1));
    waitbar(perc/100,h,sprintf('Progress:\n %0.2f%%',perc))
    comb_cell{i} = {c b};
end
close(h)
output = comb_cell;
clearvars  X x y data_dens comb_cell


% %CHANGE THIS?
% no_dens_ind = density_table.ind_nodens{g};
% temp = saved_data(no_dens_ind,:);
% temp = temp(:,{'density','sio2','tio2','al2o3','feo_tot','mgo','cao','k2o','na2o','p2o5','loi',...
%     'MALI','ASI','maficity'});

temp = data(:,{'density','sio2','tio2','al2o3','feo_tot','mgo','cao','k2o','na2o','p2o5','loi',...
    'MALI','ASI','maficity'});

temp.density_matt(:,1) = 0;
h = waitbar(0,'Initializing waitbar...');

extra_counter = 0;

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
                    %Use indices to do all of them at once so its faster?
                    ind_b = output{sum_nan}{1,1}(j,:);
                    b = output{sum_nan}{1,2}(j,:);
                    
                    %NEED AN OPTION HERE: IF THE MODEL IS EMPTY, USE THE
                    %'NEXT BEST' OPTION MODEL - check the one smaller
                    %length cell array for part matching set and use that
                    %if the b is not 0. Keep going down until sufficient.
                    
                    if ~any(b)
                        for k = 1:sum_nan-2
                            dens_avg = [];
                            flag = 0;
                            for m = 1:size(output{sum_nan-k}{1,1},1)
                                if length(intersect(find(~isnan(temp{i,2:end-2})),output{sum_nan-k}{1,1}(m,:)))==sum_nan-k
                                    ind_b = output{sum_nan-k}{1,1}(m,:);
                                    b = output{sum_nan-k}{1,2}(m,:);
                                    dens_avg(end+1) = b(1);
                                    end_ind = length(dens_avg);
                                    for l = 2:length(b)
                                        dens_avg(end_ind) = dens_avg(end_ind) + b(l)*temp{i,ind_b(l-1)+1};
                                    end
                                end
                            end

                            if any(dens_avg)
                                temp.density_matt(i,1) = dens_avg;
                                extra_counter = extra_counter+1;
                                break;
                            else
                                continue;
                            end
                        end
                        
                        
                        
                        
                    else
                    
                        temp.density_matt(i,1) = b(1);

                        for k = 2:length(b)
                        temp.density_matt(i,1) = temp.density_matt(i,1) + b(k)*temp{i,ind_b(k-1)+1};
                        end
                        break;
                    
                    end
                    
                    
                end
            end
        end
    end
    perc = 100*(i/(size(temp,1)));
    waitbar(perc/100,h,sprintf('Progress \n %0.2f%%',perc))
end



extra_counter



close(h)
temp.density_matt(temp.density_matt==0) = NaN;

%Change this
% data_save.density_matt(no_dens_ind,:) = temp.density_matt;
data_save.density_matt(data_save_ind,:) = temp.density_matt;

end



data_save.density_matt(data_save.density_matt==0) = NaN;

warning('on','all')


ind = ~isnan(saved_data.density);
derrick_data = vpest(saved_data(ind,:));
derrick_data = densest(derrick_data);
% 
% derrick_data = vpest(saved_data);
% derrick_data = densest(derrick_data);




figure()
plot(derrick_data.density,derrick_data.density_model,'.r')
hold on
plot([2000:200:4000],[2000:200:4000],'-g')
hold off
title('Derricks')
ylim([2400 3600])
xlim([2400 3600])
s = int2str(size(derrick_data.density_model(~isnan(derrick_data.density_model) & derrick_data.density_model>2000,:),1));
legend(horzcat('Number of total estimates: ',s))

figure()
plot(data_save.density,data_save.density_matt,'.b')
hold on
plot([2000:200:4000],[2000:200:4000],'-g')
hold off
title('Matts')
ylim([2400 3600])
xlim([2400 3600])
s = int2str(size(data_save.density_matt(~isnan(data_save.density_matt) & data_save.density_matt>2000,:),1));
legend(horzcat('Number of total estimates: ',s))

warning('on','all')


%output = 1;

return