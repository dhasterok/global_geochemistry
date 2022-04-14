%--------------------------------------------------------------------------
%                       Age hp plots w/ error bars
%--------------------------------------------------------------------------

%Takes in a table T with columns for age,age_min,age_max and hp

function [stand_avg] = age_hp_error(data,name_data)
fprintf('Calculating HP through time...\n')
vars = {'sample_id','age_min','age','age_max','heat_production',...
        'sio2','feo_tot','cao','mgo','k2o'};
%Remove rows with no heat_production estimate
    rows = ~isnan(data.heat_production);
    data = data(rows,vars);
    
%Correct for age_min/age_max/age missing [Add this to postgresql]
    for i = 1:size(data,1)
        if isnan(data.age(i))
            if isnan(data.age_min(i)) && ~isnan(data.age_max(i))
                data.age(i) = data.age_max(i);
            elseif ~isnan(data.age_min(i)) && isnan(data.age_max(i))
                data.age(i) = data.age_min(i);
            elseif ~isnan(data.age_min(i)) && ~isnan(data.age_max(i))
                data.age(i) = (data.age_max(i)+data.age_min(i))/2;
            else
                %Nothing
            end
        end
    end

%Remove negative age values [Add this to postgresql/a more advanced check
%that incorporates maybe just neg these values?
    rows = data.age>0;
    data = data(rows,vars);
    

%Remove NaN age, sio2, feo_tot, cao, mgo, k2o
    rows = ~isnan(data.age);
    data = data(rows,vars);
    rows = ~isnan(data.sio2);
    data = data(rows,vars);
    rows = ~isnan(data.feo_tot);
    data = data(rows,vars);
    rows = ~isnan(data.cao);
    data = data(rows,vars);
    rows = ~isnan(data.mgo);
    data = data(rows,vars);
    rows = ~isnan(data.k2o);
    data = data(rows,vars);
    
%Remove negative age values (for now)
    rows = data.age>0;
    data = data(rows,vars);

age_range = [(3800:-200:0)' (4000:-200:200)'];

stand_avg = zeros(length(age_range),6);

temp_array = cell(1,1);
for i = 1:length(age_range)
    count = 1;
    temp = [];
    for j = 1:length(data.age)
        if data.age(j)<=age_range(i,2) && data.age(j)>=age_range(i,1)
            temp(count,1) = data.heat_production(j);
            count = count+1;
        end
    end
    stand_avg(i,1) = (age_range(i,2)+age_range(i,1))./2;
    stand_avg(i,3) = length(temp);
    if isempty(temp)
        stand_avg(i,2) = NaN;
    else
        stand_avg(i,2) = sum(temp)./length(temp);
    end
    
    quant = quantile(temp,[0.25 0.50 0.75]);
    stand_avg(i,4)=quant(1,1);
    stand_avg(i,5)=quant(1,2);
    stand_avg(i,6)=quant(1,3);
    temp_array(i,:) = {temp};
%     figure()
%     h1 = histogram(temp,[0:0.5:20]);
%     hold on
%     line([quant(1,1) quant(1,1)],[0 max(h1.Values)],'Color','r');
%     line([quant(1,2) quant(1,2)],[0 max(h1.Values)],'Color','r');
%     line([quant(1,3) quant(1,3)],[0 max(h1.Values)],'Color','r');
%     hold off
%     title(stand_avg(i,1))
%     xlim([0 20])
end


% figure()
% errorbar(stand_avg(:,1),stand_avg(:,5),stand_avg(:,5)-stand_avg(:,4),stand_avg(:,6)-stand_avg(:,5))



alltemp = [];
allsizes = [];
for i = 1:length(age_range)
    alltemp = vertcat(alltemp,temp_array{i,:});
    allsizes = vertcat(allsizes,stand_avg(i,1)*ones(size(temp_array{i,:})));
end
figure('position',[300,300,1000,500])
h = boxplot(alltemp,allsizes,'Symbol','')
ylim([0,12])
xlim([0.5 length(age_range)+0.5])
s = sprintf('Age vs. HP: %s', name_data);
title(s)
xlabel('Age (Mya)') % x-axis label
ylabel('Heat Production (uW m-3)') % y-axis label

for i = 1:length(age_range)
    text(i,11,num2str(floor(stand_avg(length(age_range)-i+1,3)),'%d'),'HorizontalAlignment','center','VerticalAlignment','bottom')
end
text(length(age_range)+1.5,11,'No. of data','HorizontalAlignment','center','VerticalAlignment','bottom')
% for i=1:numel(stand_avg(:,6))
%     text(stand_avg(i,3),stand_avg(i,6),num2str(stand_avg(i,3),'%d'),...
%                'HorizontalAlignment','center',...
%                'VerticalAlignment','bottom')
% end

end