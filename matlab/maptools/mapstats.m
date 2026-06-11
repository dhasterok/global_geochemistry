function rt = mapstats(data,field,varargin)
% MAPSTATS - creates roughly equal area map of data average (median)
%
%   mapstats(data,field) where field is the field to be
%
%   mapstats(data,field,scale) where scale is 'log' or 'linear'

scale = 'linear';
dl = 2;
opt = 1;
if nargin > 2
    while opt <= nargin - 2
        switch lower(varargin{opt})
            case 'scale'
                scale = varargin{opt + 1};
                opt = opt + 2;
            case 'resolution'
                dl = varargin{opt + 1};
                opt = opt + 2;
            otherwise
                error('Unknown option.');
        end
    end
end

ind = data{:,field} > 0;
latitude = data.latitude(ind);
longitude = data.longitude(ind);
if strcmp(scale,'linear');
    d = data{ind,field}; 
elseif strcmp(scale,'log');
    d = log10(data{ind,field});
else
    error('Scale must be log or linear.');
end

rt = rtstats(dl,dl,longitude,latitude,d);
ind = ~isnan(rt.mean);
sum(rt.n>0)
global_mean = mean(rt.mean(ind))
global_std = std(rt.mean(ind))
%map = palette([0 0 0.7625],[1 0.8431 0],[0.7 0 0],64);
%map = rwb;
load('rwb_gauss.mat');


%figure;
%subplot(221);
%rtplot(rt.vert(ind,:),rt.mean(ind));
%title('(a) Median Heat Production [\muW m^{-3}]');
%caxis(log10([0.1 10.^(2*global_mean + 1)]));
%cb = colorbar;
%cb.TickDirection = 'out';
%tics = [0.1 0.3 1 10^global_mean 3 10 30];
%ticlbl = ['0.1' '0.3' '1' 'GM' '3' '10' '30'];
%cb.Ticks = log10(tics);
%cb.TickLabels = tics;
%colormap(map);
figure;
subplot(211);
norm_mean = (rt.mean(ind) - global_mean)/global_std;
rtplot(rt.vert(ind,:),norm_mean);
cb.TickDirection = 'out';

%if strcmp(field,'heat_production')
    caxis([-3 3]);
    if strcmp(scale,'linear')
        [tics,ord] = sort([-3:0.5:3]);

        tlbl = {'-4\sigma','','-3\sigma','','-2\sigma','','-1\sigma','','0\sigma, GM', ...
            '','1\sigma','','2\sigma','','3\sigma','','4\sigma'};
    else
        [tics,ord] = sort([-3:3 (log10([0.1 0.3 1 3 10 30])-global_mean)/global_std]);

        tlbl = {'-4\sigma','-3\sigma','-2\sigma','-1\sigma','0\sigma, GM', ...
            '1\sigma','2\sigma','3\sigma','4\sigma', ...
            '0.1','0.3','1','3','10','30'};
    end    
    cb.TickLabels = tlbl(ord);
    
    title('(a) Heat production [\muW m^-3] and Normalized difference from global mean (in log-space).');
%else
%    caxis([-1 1]);
%    [tics,ord] = sort([-2:0.5:2]);

    %tlbl = {'-4\sigma','','-3\sigma','','-2\sigma','','-1\sigma','','0\sigma, GM', ...
    %    '','1\sigma','','2\sigma','','3\sigma','','4\sigma'};
    
%    title(['(a) ,',field,' and Normalized difference from global mean (in log-space).']);
%end
cb.Ticks = tics;    
colormap(map);
cb = colorbar;


subplot(212);
rtplot(rt.vert(ind,:),log10(rt.n(ind)));
caxis([0 3]);
cb = colorbar;
cb.Ticks = [0:3];
cb.TickLabels = {'1','10','100','1000'};
cb.TickDirection = 'out';
title('(b) Number of data in each ~equal area bin.');

figure;
h = histogram(rt.mean(ind));
hold on;
histogram(rt.mean(ind & rt.vert(:,1) > 110 & rt.vert(:,2) < 155 & rt.vert(:,3) > -45 & rt.vert(:,4) < -10),'BinEdges',h.BinEdges);
aus_mean = mean(rt.mean(ind & rt.vert(:,1) > 110 & rt.vert(:,2) < 155 & rt.vert(:,3) > -45 & rt.vert(:,4) < -10));
yl = get(gca,'YLim');
f = fill(global_mean + [global_std global_std -global_std -global_std],[yl fliplr(yl)],[0.5 0.5 0.5]);
f.EdgeColor = 'none';
f.FaceAlpha = 0.3;
plot([global_mean,global_mean],yl,'k-');
plot([aus_mean,aus_mean],yl,'k--');
ylim(yl);
if strcmp(scale,'log')
    hpax([floor(min(rt.mean)) ceil(max(rt.mean))],'x');
end
xlabel(field);
title('(c) Distribution of equal area means.');
golden

return