function chemplot(data,x,xl,xunit,y,yl,yunit,varargin)

if any(strcmpi(x,data.Properties.VariableNames))
    ind = find(strcmpi(x,data.Properties.VariableNames));
    X = data{:,ind};
else
    warning(['Could not find data field ',x]);
end

if length(xl) == 3
    if xl(3) == 1
        X = log10(X);
    end
end

if any(strcmpi(y,data.Properties.VariableNames))
    ind = find(strcmpi(y,data.Properties.VariableNames));
    Y = data{:,ind};
else
    warning(['Could not find data field ',y]);
end

if length(yl) == 3
    if yl(3) == 1
        Y = log10(Y);
    end
end

%ind = d > 0;
if nargin == 4
    ind = ~isnan(X) & ~isnan(Y) & varargin{1};
else
    ind = ~isnan(X) & ~isnan(Y);
end
ex = linspace(xl(1),xl(2),40);
ey = linspace(yl(1),yl(2),40);

n = hist2d(X(ind),Y(ind),ex,ey);
n(n == 0) = NaN;
n = log10(n);
%[c,h] = contourf(ec(1:end-1) + diff(ec)/2,ea(1:end-1) + diff(ea)/2,n);
imagesc(ex(1:end-1) + diff(ex)/2,ey(1:end-1) + diff(ey)/2,n);
colormap(flipud(gray));
caxis([0 3]);
axis xy;

if length(xl) == 2
    xlabel([x,' ',xunit]);
    xlim(xl);
elseif xl(3) == 0
    xlabel([x,' ',xunit]);
    xlim(xl);
else
    hpax(xl(1:2),'x');
    xlabel([x,' ',xunit]);
end


if length(yl) == 2
    ylabel([y,' ',yunit]);
    ylim(yl);
elseif yl(3) == 0
    ylabel([y,' ',yunit]);
    ylim(yl);
else
    hpax(yl(1:2),'y');
    ylabel([y,' ',yunit]);
end

return

