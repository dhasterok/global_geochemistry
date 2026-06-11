function spider_numbers(t)
% SPIDER_NUMBERS Plots number of analyses in a spider plot
%
%   spider_numbers(t) plots the number of analyses in each element for each
%   series/data subset in a spider datagram.  t must be a table generated
%   by spider.

% find column with number of data
ind = contains(t.Properties.VariableNames,'_N');

% subsets, y-axis label
subset = t.subset;
% elements, x-axis label
ellist = t.Properties.VariableNames(ind);
% convert column names to elements
for i = 1:length(ellist)
    ellist{i} = ellist{i}(1:end-6);
    ellist{i}(1) = upper(ellist{i}(1));
end

% matrix of number of data
N = t{:,ind};

% plot numbers
imagesc(log10(N));
golden

% axes labels
set(gca,'TickDir','out', ...
    'Ytick',[1:length(subset)],'YTickLabel',subset, ...
    'XTick',[1:length(ellist)],'XTickLabel',ellist);

% colorbar
cb = colorbar;
cb.TickDirection = 'out';

clim = cb.Limits;
caxis([0 ceil(clim(2))]);
cb.Ticks = [0:ceil(clim(2))];
cb.TickLabels = 10.^cb.Ticks;
cb.Label.String = 'No. Analyses';

return