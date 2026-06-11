function [output density_matt] = density_est(data)
% lastwarn('')
warning('off','all')
data.density(data.density > 0 & data.density < 10,:) = data.density(data.density > 0 & data.density < 10,:)*1000;

data = data(data.sio2 > 0 & data.sio2 < 100,:);
data = data((data.al2o3 > 0 & data.al2o3 < 100) | isnan(data.al2o3),:);
data = data((data.feo_tot > 0 & data.feo_tot < 100) | isnan(data.feo_tot),:);
data = data((data.mgo > 0 & data.mgo < 100) | isnan(data.mgo),:);
data = data((data.cao > 0 & data.cao < 100) | isnan(data.cao),:);
data = data((data.na2o > 0 & data.na2o < 100) | isnan(data.na2o),:);
data = data((data.k2o > 0 & data.k2o < 100) | isnan(data.k2o),:);

y = data.density;
x{:,1} = data.sio2;
x{:,2} = data.al2o3;
x{:,3} = data.feo_tot;
x{:,4} = data.mgo;
x{:,5} = data.cao;
x{:,6} = data.na2o;
x{:,7} = data.k2o;
x{:,8} = data.Fe_number;
x{:,9} = data.MALI;
x{:,10} = data.ASI;
x{:,11} = data.maficity;

n_var = size(x,2);
comb_cell = cell(1,n_var-1);

%X = [ones(size(x1)) x1 x2 x3 x4 x5 x6 x7];
%b = regress(y,X);
%loop through: estimate = b, b(1) + b(2)*x1 + ... b(n)*xn

h = waitbar(0,'Initializing waitbar...');

data_dens = data(data.density>0,:);

for i = 1:n_var-1
    %Combinations must include sio2
    c = [ones(size(combnk(2:n_var,i),1),1) combnk(2:n_var,i)];
    b = [];
    for j = 1:size(c,1)
        ind = [];
        y = data_dens.density;
        x{:,1} = data_dens.sio2;
        x{:,2} = data_dens.al2o3;
        x{:,3} = data_dens.feo_tot;
        x{:,4} = data_dens.mgo;
        x{:,5} = data_dens.cao;
        x{:,6} = data_dens.na2o;
        x{:,7} = data_dens.k2o;
        x{:,8} = data_dens.Fe_number;
        x{:,9} = data_dens.MALI;
        x{:,10} = data_dens.ASI;
        x{:,11} = data_dens.maficity;
        
        for k = 1:size(c,2)
            ind{k} = find(~isnan(x{:,c(j,k)}));
        end
        
        if isempty(ind)
            continue;
        end
        
        com_ind = ind{1};
        for k = 2:size(c,2)
            com_ind = intersect(com_ind,ind{k});
        end

        ind = com_ind;
        
        if isempty(ind)
            continue;
        end
        
        y = data_dens.density(ind,:);
        x{:,1} = data_dens.sio2(ind,:);
        x{:,2} = data_dens.al2o3(ind,:);
        x{:,3} = data_dens.feo_tot(ind,:);
        x{:,4} = data_dens.mgo(ind,:);
        x{:,5} = data_dens.cao(ind,:);
        x{:,6} = data_dens.na2o(ind,:);
        x{:,7} = data_dens.k2o(ind,:);
        x{:,8} = data_dens.Fe_number(ind,:);
        x{:,9} = data_dens.MALI(ind,:);
        x{:,10} = data_dens.ASI(ind,:);
        x{:,11} = data_dens.maficity(ind,:);
        
        X = [ones(size(x{:,1}))];
        for k = 1:size(c,2)
            X = [X x{:,c(j,k)}];
        end
        
        b(j,:) = regress(y,X)';

        perc = 100*(i/(n_var-1)) + 10*(j/size(c,1));
        waitbar(perc/100,h,sprintf('Progress \n %0.2f%%',perc))
%         if ~isempty(lastwarn)
%             c(j,:)
%             y
%             [x{:,1} x{:,5} x{:,6} x{:,7} x{:,9}]
%             b(j,:)
%             return
%         end
%         if isequal(c(j,:),[1 2 4 5 6 7 9])
%             c(j,:)
%             b(j,:)
%             return
%         end
    end
    comb_cell{i} = {c b};
end

close(h)

output = comb_cell;

clearvars  X x y data_dens comb_cell

temp = data(:,{'density','sio2','al2o3','feo_tot','mgo','cao','na2o','k2o',...
    'Fe_number','MALI','ASI','maficity','density_model'});

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

figure()
plot(temp.density,temp.density_model,'.r')
hold on
plot([2000:200:4000],[2000:200:4000],'-g')
hold off
title('Derricks')
ylim([2400 3600])
xlim([2400 3600])
s = int2str(size(temp.density_model(~isnan(temp.density_model),:),1));
legend(horzcat('Number of total estimates: ',s))

figure()
plot(temp.density,temp.density_matt,'.b')
hold on
plot([2000:200:4000],[2000:200:4000],'-g')
hold off
title('Matts')
ylim([2400 3600])
xlim([2400 3600])
s = int2str(size(temp.density_matt(~isnan(temp.density_matt),:),1));
legend(horzcat('Number of total estimates: ',s))

density_matt = temp.density_matt;

warning('on','all')

return

% data_1234567891011 = data(~isnan(data.sio2) & ~isnan(data.al2o3)...
%      & ~isnan(data.feo_tot) & ~isnan(data.mgo) & ~isnan(data.cao)...
%       & ~isnan(data.na2o) & ~isnan(data.k2o) & ~isnan(data.Fe_number)...
%        & ~isnan(data.MALI) & ~isnan(data.ASI) & ~isnan(data.maficity),:);