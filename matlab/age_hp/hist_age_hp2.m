%Data - sio2, feo, mgo, cao, k2o
data2 = [randn(1,100); randn(1,100); randn(1,100); randn(1,100); randn(1,100);...
    randn(1,100);randn(1,100);randn(1,100);randn(1,100);randn(1,100);...
    randn(1,100);randn(1,100);randn(1,100);randn(1,100);randn(1,100);...
    randn(1,100);randn(1,100);randn(1,100);randn(1,100);randn(1,100)];

for i = 1:size(data2,1)
    h = histogram(data2(i,:),-3:0.5:3);
    prop = properties(h);
    for j = 1:length(prop)
        newh(i).(prop{j}) = h.(prop{j});
    end
    newh(i).BarEdges = (newh(i).BinEdges(:,1:end-1)+newh(i).BinEdges(:,2:end)) ./ 2;
end

close all
figure()
subplot(2,1,1)

xlocations = [1 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21];
%xlocations = 1:size(data,1);
scatter(xlocations,1:size(data2,1),'o');
xlim([0.5 max(xlocations)+0.5])

%Get position of a point
h = subplot(2,1,1);
pos = get(h,'position');
xlim_var = get(h,'xlim');
ylim_var = get(h,'ylim');
for i = 1:size(data2,1)
    x_in_pixels(i) = pos(1) + pos(3) * (xlocations(i)-xlim_var(1))/(xlim_var(2)-xlim_var(1));
end

min_diff_x = min(x_in_pixels(:,2:end)-x_in_pixels(:,1:end-1));

for i = 1:(size(newh,2))
    ax_hists(i) = axes('Position',[x_in_pixels(i)-min_diff_x/4 0.1 min_diff_x/2 0.4]);
    %Something
    hold on
    barh(newh(i).BarEdges,newh(i).Values,1)
    barh(newh(i).BarEdges,-newh(i).Values,1)
    hold off

    %USE AXES INSTEAD OF SUBPLOT IF SETTING POSITIONS
    %has a lower left corner at the point (0.1 0.1) with a width and height of 0.7
%     figure
%     ax1 = axes('Position',[0.1 0.1 0.7 0.7]);
%     ax2 = axes('Position',[0.65 0.65 0.28 0.28]);
    
    
    
    
    %Make these 25-85 or something for SiO2, and x lim (weights) need to be
    %different, only want to show proportions. Dont put any restrict on
    %this then?
end