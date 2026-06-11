function feox(data,varargin)

rtype = '';
if nargin == 2
    rtype = varargin{1};
end

if ~isempty(rtype)
    ind = strcmp(rtype,data.rock_type);
else
    ind = logical(ones(size(data.FeO)));
end
ind = ind & data.FeO > 0 & data.Fe2O3 > 0 & data.heat_production > 0;
nfe2 = data.FeO(ind)/71.8444;
nfe3 = 2*data.Fe2O3(ind)/159.6882;

oxstate = nfe3./(nfe3 + nfe2) + 2;

figure;
ec = [0:0.1:15];
ea = [-3:0.1:2];
subplot(221);
n = log10(hist2d(data.FeO(ind),log10(data.heat_production(ind)),ec,ea));
imagesc(ec(1:end-1) + diff(ec)/2,ea(1:end-1) + diff(ea)/2,n);
colormap(flipud(gray));
hpax([-3 2]);
xlabel('FeO (wt.%)');
axis xy;

subplot(223);
ec = [0:0.1:8];
ea = [-3:0.1:2];
n = log10(hist2d(data.Fe2O3(ind),log10(data.heat_production(ind)),ec,ea));
imagesc(ec(1:end-1) + diff(ec)/2,ea(1:end-1) + diff(ea)/2,n);
colormap(flipud(gray));
hpax([-3 2]);
xlabel('Fe2O3 (wt.%)');
axis xy;

subplot(222);
ec = linspace(2,3,40);
ea = [-3:0.1:2];
n = log10(hist2d(oxstate,log10(data.heat_production(ind)),ec,ea));
imagesc(ec(1:end-1) + diff(ec)/2,ea(1:end-1) + diff(ea)/2,n);
colormap(flipud(gray));
%plot(oxstate,data.heat_production(ind),'.');
hpax([-3 2]);
xlabel('Iron oxidation state, Fe^{n+}');
axis xy;

subplot(224);
plot(data.SIO2(ind),oxstate,'.');
ec = [30:0.2:95];
ea = [2:0.02:3];
n = log10(hist2d(data.SIO2(ind),oxstate,ec,ea));
imagesc(ec(1:end-1) + diff(ec)/2,ea(1:end-1) + diff(ea)/2,n);
colormap(flipud(gray));
xlim([35 90]);
xlabel('SiO2 [wt.%]');
ylabel('Iron oxidation state, Fe^{n+}');
axis xy;

return
