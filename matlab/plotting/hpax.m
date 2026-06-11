function hpax(lim,varargin)
% hpax(lim) where lim are the log10 values of the axes limits
%
% hpax(lim,x_or_y) where x_or_y is either 'x' or 'y' for axis (y is assumed)

if nargin == 1
    x_or_y = 'y';
elseif nargin == 2
    x_or_y = varargin{1};
else
    help hpax;
end

mt = log10([1:9]);
tmp = {'','','','','','','',''};
atic = [];
aticlbl = {};
for i = lim(1):lim(2)
    atic = [atic,i + mt];
    aticlbl = [aticlbl,{num2str(10^i)},tmp];
end
if strcmp(x_or_y,'x')
    xlabel('Heat Production [µW m^{-3}]');
    %set(gca,'XTick',[lim(1):lim(2)],'XTickLabel',10.^[lim(1):lim(2)],'Box','on');
    set(gca,'XTick',atic,'XTickLabel',aticlbl,'Box','on');
    xlim(lim);
elseif strcmp(x_or_y,'y')
    ylabel('Heat Production [µW m^{-3}]');
    %set(gca,'YTick',[lim(1):lim(2)],'YTickLabel',10.^[lim(1):lim(2)],'Box','on');
    set(gca,'YTick',atic,'YTickLabel',aticlbl,'Box','on');
    ylim(lim);
else
    warning('Ignoring hpax command, incorrect axes argument');
end

return
