% Data
X=[data.sio2, data.k2o + data.na2o,data.heat_production];

% Compute the histogram
edge = linspace(0,140,200); % change here to your need;
[count, ~, ~, loc] = histcn(X,edge,edge,edge);

% Gaussian smoothing the histogram
kernel = exp(-linspace(-2,2,11).^2);
K = 1;
for k=1:3
    K = K(:)*kernel;
end
K = reshape(K,length(kernel)+[0 0 0 ]);
K = K/sum(K(:));
count = convn(count, K,'same');

% Density
density= zeros(size(X,1),1);
valid = all(loc,2);
loc = loc(valid,:);
density(valid) = count(sub2ind(size(count),loc(:,1),loc(:,2),loc(:,3)));

density = density;

% Plot
scatter3(X(:,1),X(:,2),X(:,3),5,density)
colormap(hot)


set(gca,'XLim',[30 90],'YLim',[0 15],'ZLim',[0 10])
xlabel('sio2')
ylabel('TA')
zlabel('hp')