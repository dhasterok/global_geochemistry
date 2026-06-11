function [coeff, score, latent, Tsquared, explained, el,X,avg] = hp_pca_ing(data,elements)

% X - data vector for PCA [each sample, each element/oxide]
% since not all oxides and elements are measured for every sample, we need
% a set of oxides and elements that are found in the majority of analyses

%cull data to only get rocks in the lower mantle
data = data(rockgroup(data,'all plutonic'),:);

for i = 1:length(elements)
    n(i) = sum(~isnan(data{:,elements{i}}) & data{:,elements{i}} > 0);
end

figure;
barh(n/height(data));
for i = 1:length(elements)
    if ~isempty(strfind(elements{i},'_ppm'))
        ellab{i} = elements{i}(1:end-4);
    else
        ellab{i} = elements{i};
    end
end
set(gca,'YTick',[1:length(n)],'YTickLabel',ellab);
xlabel('Number of positive values');
ylim([0 length(n)+1]);
axis ij;

thresh = 0.25;

el_use = n/height(data) > thresh;

fprintf('Using the following for PCA based on a threshold of %f\n',thresh);
for i = 1:length(el_use)
    if el_use(i)
        fprintf(' %s\n',ellab{i});
    end
end

use = ones([height(data),1]);
for i = 1:length(el_use)
    if el_use(i)
        use = use & data{:,elements{i}} > 0;
    end
end

el = ellab(el_use);

% X - data vector for PCA [each sample, each element/oxide]
X = zeros([sum(use) sum(el_use)]);
c = 1;
for i = 1:length(el_use)
    if el_use(i)
        if isempty(strfind(elements{i},'ppm')) ...
                & ~strcmp('k2o',elements{i})
            X(:,c) = data{use,elements{i}};
        else
            X(:,c) = log(data{use,elements{i}});
        end
        c = c + 1;
    end
end

%Check for bad values
X(X<0) = 0;
%ind=x>0&y>0
%data=data(ind)
avg=mean(X)
out_put1=cov(X)

% coeff - principle component coefficient
% score - principle component score
% latent - principle component variances
% Tsquared - Hotelling's T-squared statistic
% explained - percentage of total variance explained by each principle
%       component
[coeff, score, latent, Tsquared, explained] = pca(X);

coeff
latent'
explained'

figure;
subplot(121);
plot(log10(explained),'o');
xlim([0 length(el)+1]);
set(gca,'XTick',[1:length(el)]);
xlabel('Principle Component');
hpax([-2 2]);
ylabel('Percentage of total variance');
axis square;

subplot(122);
imagesc(coeff);
xlabel('Principle component coefficient');
set(gca,'YTick',[1:length(el)],'YTickLabel',el);
set(gca,'XTick',[1:length(el)])
map = load('bwr.map');
colormap(map);
cbar;
caxis([-1 1]);
axis square;

pcplot(data,score,use,'rock_origin','plutonic');

edgeStart = -40;
dx = 1;
N = 80;
edges = edgeStart + (0:N-1)*dx;
figure
for i=1:4
    subplot(4,1,i)
    hist(score(:,i),edges)
    xlim([-40,40]);
    num = num2str(i);
    T = ['Princible Axis ' num];
    title(T)
    ylabel('# data points')
    xlabel('score')
end
figure
for i=5:7
    subplot(3,1,(i-4))
    hist(score(:,i),edges)
    xlim([-40,40]);
    num = num2str(i);
    T = ['Princible Axis ' num];
    title(T)
    ylabel('# data points')
    xlabel('score')
end
%dissimilarities = pdist(zscore(X),'cityblock');
%[Y,stress] =... 
%mdscale(dissimilarities,2,'criterion','metricstress');
%plot(Y(:,1),Y(:,2),'o','LineWidth',2);
%whos


return

function pcplot(data,score,use,property,ttxt);

figure;
rtype = lower(data{use,property});

%ax3 = [-40 60];
%ax2 = [-40 60];
%ax1 = [-100 50];
ax3 = [-30 40];
ax2 = [-30 40];
ax1 = [-70 50];
x1 = linspace(ax1(1),ax1(2),50);
x2 = linspace(ax2(1),ax2(2),50);
x3 = linspace(ax3(1),ax3(2),50);


subplot(222);
ind = strcmp(lower(ttxt),rtype);
%p = plot(score(ind,1),score(ind,2),'o');
%set(p,'MarkerFaceColor',get(p,'Color'));
n = hist3([score(ind,1),score(ind,2)],'Edges',{x1 x2});
imagesc(x1,x2,n');%log10(n'));
colormap(flipud(gray));
%caxis([-0.1 3]);
xlabel('PC1');
ylabel('PC2');
%title(ttxt);
axis([ax1 ax2]);
axis square;
set(gca,'TickDir','out');

subplot(224);
%p = plot(score(ind,1),score(ind,3),'o');
%set(p,'MarkerFaceColor',get(p,'Color'));
n = hist3([score(ind,1),score(ind,3)],'Edges',{x1 x3});
imagesc(x1,x3,n');%log10(n'));
%caxis([-0.1 3]);
xlabel('PC1');
ylabel('PC3');
axis([ax1 ax3]);
axis square;
set(gca,'TickDir','out');

subplot(221);
%p = plot(score(ind,3),score(ind,2),'o');
%set(p,'MarkerFaceColor',get(p,'Color'));
n = hist3([score(ind,3),score(ind,2)],'Edges',{x3 x2});
imagesc(x3,x2,n');%log10(n'));
%caxis([-0.1 3]);
xlabel('PC3');
ylabel('PC2');
axis([ax3 ax2]);
axis square;
set(gca,'TickDir','out');

return