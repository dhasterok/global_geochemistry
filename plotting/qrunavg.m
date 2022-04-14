function [Xm,Y,N] = qrunavg(x,y,dx,varargin)

q = [0.025 0.25 0.5 0.75 0.975];
%X = [min(x)+0.5*dx:dx:max(x)-0.5*dx]';
X = [0+0.5*dx:dx:max(x)-0.5*dx]';
Xm = midpt(X);

Y = zeros([length(X)-1 length(q)]);
N = zeros([length(X)-1 1]);
for i = 1:length(X)-1
    ind = X(i) <= x & x < X(i+1);
    N(i,1) = sum(ind);
    Y(i,:) = quantile(y(ind),q);
end

%figure;
subplot(3,1,1:2); hold on;
% 95% bounds
p{1} = plot(Xm,Y(:,[1,5]),'--','LineWidth',0.25,'Color',[0.7 0.7 0.7]);
% 50% bounds
if nargin == 4
    C = varargin{1}+0.2;
    C(C>1) = 1;
    D = C+0.4;
    D(D>1) = 1;
    fill([Xm; flipud(Xm)],  [Y(:,2); flipud(Y(:,4))], ...
        D,'EdgeColor','None','FaceAlpha',0.3);
else
    fill([Xm; flipud(Xm)],  [Y(:,2); flipud(Y(:,4))], ...
        [0.5 0.5 0.5],'EdgeColor','None','FaceAlpha',0.3);
end
% 2-SE
if nargin == 4
    fill([Xm; flipud(Xm)], ...
        [Y(:,3); flipud(Y(:,3))] + [(Y(:,1)-Y(:,3)); flipud(Y(:,5)-Y(:,3))]./sqrt([N; flipud(N)]), ...
        C,'EdgeColor','None','FaceAlpha',0.3);
else
    fill([Xm; flipud(Xm)], ...
        [Y(:,3); flipud(Y(:,3))] + [(Y(:,1)-Y(:,3)); flipud(Y(:,5)-Y(:,3))]./sqrt([N; flipud(N)]), ...
        [1 0 0],'EdgeColor','None','FaceAlpha',0.3);    
end
% median
p{2} = plot(Xm,Y(:,3),'r-','LineWidth',1.5);
if nargin == 4
    for i = 1:length(p)
        set(p{i},'Color',varargin{1});
    end
end
set(gca,'Box','on');
golden;

pba = get(gca,'PlotBoxAspectRatio');

subplot(3,1,3); hold on;

if nargin == 4
    histogram(x,'DisplayStyle','stairs');
    %set(gca,'Color',varargin{1});
else
    histogram(x,'DisplayStyle','stairs');
    %set(gca,'Color',[1 0 0]);
end

set(gca,'PlotBoxAspectRatio',[2.38*pba(1) pba(2:3)],'Box','on');

return