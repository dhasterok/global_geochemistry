function density_check(data,modelfield)

% observed density
density_obs = data.density;
% make sure in kg m^-3
ind = density_obs < 10;
density_obs(ind) = density_obs(ind)*1000;

% modeled density
density_model = data{:,{modelfield}};
ind = ~strcmp('carbonatite',data.rock_name) & ~strcmp('',data.rock_name);

% range of plot
ax = [2400 3600];

% compute 2D histogram of observed vs. modeled density
er = [ax(1):15:ax(2)]; % edges
n = hist2d(density_obs(ind),density_model(ind),er,er);

%plot(density(ind),density_model(ind),'.');
% plot histogram as image
imagesc(er(1:end-1) + diff(er)/2,er(1:end-1) + diff(er)/2,log10(n));

colormap(flipud(gray));
axis xy;
% logarithmic colorscale
caxis([-0.1 2]);
colorbar;

% compute statistics for model misfit
r = (density_obs(ind) - density_model(ind));
r = r(~isnan(r));
fprintf('Average difference between models: %0f\n',mean(r));
sd = std(r);

% plot 1:1 line
hold on;
plot(ax,ax,'-');
text(ax(1)+100,ax(2)-100,['\sigma = ',num2str(sd)]);
axis([ax ax]);
axis square;
xlabel('Observed Density (kg m^{-3})');
ylabel('Estimated Density (kg m^{-3})');
title(modelfield);

return