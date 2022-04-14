function chemhpplot(data,el,xl,varargin)

if any(strcmpi(el,data.Properties.VariableNames))
    eind = find(strcmpi(el,data.Properties.VariableNames));
    d = data{:,eind};
else
    warning(['Could not find data field ',el]);
end
%ind = d > 0;
if nargin == 4
    ind = ~isnan(d) & varargin{1};
else
    ind = ~isnan(d);
end
ec = linspace(xl(1),xl(2),40);
ea = [-3:0.1:2];

n = log10(hist2d(d(ind),log10(data.heat_production(ind)),ec,ea));
%[c,h] = contourf(ec(1:end-1) + diff(ec)/2,ea(1:end-1) + diff(ea)/2,n);
imagesc(ec(1:end-1) + diff(ec)/2,ea(1:end-1) + diff(ea)/2,n);
colormap(flipud(gray));
caxis([0 3]);
axis xy;
xlabel([el,' (wt. %)']);
xlim(xl);
hpax([-3 2]);

return
